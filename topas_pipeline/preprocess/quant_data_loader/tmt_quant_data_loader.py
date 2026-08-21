import logging
import numpy as np
import pandas as pd
from typing import List, Dict

from topas_pipeline.preprocess.file_processor.result_file_loader_factory import (
    ResultFileLoaderFactory,
)
from topas_pipeline.preprocess.preprocess_tools import (
    get_save_correction_factors_function,
    get_save_debug_df_function,
    save_qc_info,
)
from topas_pipeline.preprocess.quant_data_loader.quant_data_loader import (
    BaseQuantDataLoader,
)
from topas_pipeline import utils
from topas_pipeline.sample_annotation import load_sample_annotation
from topas_pipeline.preprocess.file_processor.constants import MQ_INPUT_COLUMNS_TMT

logger = logging.getLogger(__name__)

class TMTQuantDataLoader(BaseQuantDataLoader):
    def __init__(
        self,
        results_file_list: List[str],
        quant_file_format: str,
        sample_annotation_df: pd.DataFrame,
        normalize_to_reference: bool,
        debug: bool,
        metadata: Dict,
    ):
        self.quant_file_format = quant_file_format
        self.sample_annotation_df = sample_annotation_df
        self.normalize_to_reference = normalize_to_reference
        self.debug = debug
        self.metadata = metadata
        self.results_folder = self.metadata["results_folder"]
        self.data_type = self.metadata["data_type"]
        self.quant_strategy = self.metadata["quant_strategy"]

        # Get Result File loader based on quant file format used (Maxquant, Ionquant, Diann)
        result_file_loader_class = ResultFileLoaderFactory.get_loader(
            "_".join([self.quant_strategy, self.quant_file_format])
        )
        self.result_file_loader = result_file_loader_class(
            results_file_list, sample_annotation_df=self.sample_annotation_df
        )

    def load_and_normalize(self) -> pd.DataFrame:
        """
        Load and normalize the data using the appropriate ResultFileLoader
        and return the DataFrame with normalized values
        """
        # Load the data using the load method
        df = self.load()
        save_qc_info(df, self.sample_annotation_df, self.results_folder, self.data_type)

        # Normalize the data using the normalize method
        ref_channels_df = self.sample_annotation_df.loc[
            self.sample_annotation_df["is_reference"], :
        ]
        df = self.normalize(df, ref_channels_df, self.normalize_to_reference)
        return df

    def load(self) -> pd.DataFrame:
        """
        Load the result files using the appropriate ResultFileLoader
        and return a concatenated DataFrame
        """
        # Load the result file using the load function of the ResultFileLoader
        all_batches = self.result_file_loader.load(usecols=MQ_INPUT_COLUMNS_TMT)
        df = pd.concat(all_batches, ignore_index=True)
        return df

    def normalize(
        self,
        df: pd.DataFrame,
        ref_channels_df: pd.DataFrame,
        normalize_to_reference: bool,
    ) -> pd.DataFrame:
        """
        Normalize the intensities using median centering and MS1 normalization,
        and return the DataFrame
        """
        # Get results folder and data type from metadata to save intermediate results
        results_folder = self.metadata["results_folder"]
        data_type = self.metadata["data_type"]
        save_debug_df = get_save_debug_df_function(
            self.debug, results_folder, data_type
        )
        save_correction_factors = get_save_correction_factors_function(
            results_folder, data_type
        )
        save_debug_df(df, "_complete_raw")

        # 1) Medien centering within batch
        df, correction_factors = self.median_centering_within_batch(df)
        save_debug_df(df, "_after_1st_median")
        save_correction_factors(correction_factors, "_in_batch_correction_factors")

        # 2) Medien centering MS1
        df, correction_factors = self.median_centering_ms1(df)
        save_debug_df(df, "_after_ms1_centering")
        save_correction_factors(correction_factors, "_ms1_correction_factors")

        # 3) If normalize to ref is true: scale ms1 to ref else: impute ms1 intensity
        #       requires ref_channel_df for if else in step 3
        if normalize_to_reference:
            # Scale MS1 intensities such that reference channel MS1 contribution is constant across batches
            df = self.scale_ms1_to_reference(df, ref_channels_df)
            save_debug_df(df, "_after_ms1_normalize_to_reference")
        else:
            # Impute an MS1 value for scans without MS1
            # This uses the assumption that the reference channel abundance should be constant across batches
            df = self.impute_ms1_intensity(df, ref_channels_df)
            save_debug_df(df, "_after_ms1_imputation")

        # 4) Redistribute ms1 intensity
        # Distribute MS1 intensity over the TMT channels
        df = self.redistribute_ms1_intensity(df)
        save_debug_df(df, "_after_ms1_correction")

        return df
    
    def scale_ms1_to_reference(self, df: pd.DataFrame, ref_channel_df: pd.DataFrame) -> pd.DataFrame:
        """Scale MS1 intensities such that reference channel MS1 contribution is constant across batches.

        TODO: combine this with _impute_ms1_intensity.

        For each modified peptide in each batch check how many batches have MS1 and reference channels.
        For those, compute the mean intensity of the reference channels after MS1 redistribution.
        :param df: dataframe of batch intensities (11 reporter channels, intensities are not log transformed and missing values as NaNs
                                                MS1 column, intensities are not log transformed and missing values as zeroes)
        """
        logger.info("Scale MS1 intensities using reference channels")
        df["MS1"] = df["MS1"].replace(0, np.nan)

        tmt_channels = utils.filter_for_intensity_columns(df)
        ref_channels_mean = utils.get_ref_channels_mean(df, ref_channel_df)

        df["MS1 corrected reference intensity"] = (
            ref_channels_mean / tmt_channels.sum(axis=1) * df["MS1"]
        )

        # TODO: replace this with sorted string concat to not get accidental collisions, e.g. (lot1+lot3)/2 != lot2
        batch_to_qc_lot_mapping = ref_channel_df.groupby("Batch")["QC Lot"].agg(
            "mean"
        )

        df["QC Lot"] = df["Batch"].map(batch_to_qc_lot_mapping)

        imputation_df = (
            df.groupby(["Modified sequence", "Charge", "QC Lot"])["MS1 corrected reference intensity"]
            .apply(lambda x: x.mean())
            .reset_index(name="Mean MS1 corrected reference intensity")
        )

        df = pd.merge(
            left=df, right=imputation_df, on=["Modified sequence", "Charge", "QC Lot"], how="left"
        )

        has_reference = ~np.isnan(ref_channels_mean)
        df.loc[has_reference, "MS1"] = (
            df["Mean MS1 corrected reference intensity"]
            * tmt_channels.sum(axis=1)
            / ref_channels_mean
        )
        df.loc[~has_reference, "MS1"] = np.nan

        return df
    
    @staticmethod
    def impute_ms1_intensity(df: pd.DataFrame, ref_channel_df: pd.DataFrame) -> pd.DataFrame:
        """Impute an MS1 value for scans without MS1 but with reference channel values.
        For each modified peptide in each batch if MS1 missing check how many other batches have MS1 and reference channels.
        For those, compute the mean intensity of the reference channels after MS1 redistribution.
        :param df: dataframe of batch intensities (11 reporter channels, intensities are not log transformed and missing values as NaNs
                                                MS1 column, intensities are not log transformed and missing values as zeroes)
        """
        logger.info("Imputing missing MS1 intensities using reference channels")
        df["MS1"] = df["MS1"].replace(0, np.nan)

        tmt_channels = utils.filter_for_intensity_columns(df)
        ref_channels_mean = utils.get_ref_channels_mean(df, ref_channel_df)

        df["MS1 corrected reference intensity"] = (
            ref_channels_mean / tmt_channels.sum(axis=1) * df["MS1"]
        )

        # TODO: switch to ['Modified sequence', 'Charge', 'QC Lot']
        imputation_df = (
            df.groupby("Modified sequence")["MS1 corrected reference intensity"]
            .apply(lambda x: x.mean())
            .reset_index(name="Mean MS1 corrected reference intensity")
        )

        df = pd.merge(left=df, right=imputation_df, on="Modified sequence", how="left")

        missing_ms1 = np.isnan(df["MS1"])
        has_reference = ~np.isnan(ref_channels_mean)
        df.loc[missing_ms1 & has_reference, "MS1"] = (
            df["Mean MS1 corrected reference intensity"]
            * tmt_channels.sum(axis=1)
            / ref_channels_mean
        )

        # TODO: impute MS1 for scans with no MS1 and no reference channel values by using summed channel intensities (seems only necessary for FP)

        num_ms1_before = (~missing_ms1).sum()
        num_imputed = (
            missing_ms1
            & has_reference
            & ~np.isnan(df["Mean MS1 corrected reference intensity"])
        ).sum()
        logger.info(f"Found {num_ms1_before} scans with MS1.")
        logger.info(
            f"Imputed MS1 values for {num_imputed} additional scans (+{round(num_imputed / num_ms1_before * 100, 1)}%)"
        )

        return df


    @staticmethod
    def redistribute_ms1_intensity(df: pd.DataFrame) -> pd.DataFrame:
        """Multiply MS1 with share of intensity in the channel relative to the summed intensity.
        :param df: dataframe of batch intensities (11 reporter channels, intensities are not log transformed and missing values are NaNs
                                                MS1 column, intensities are not log transformed and missing values are zeroes)
        :return:
        """
        logger.info("Distributing MS1 intensity over TMT channels")
        has_ms1 = ~np.isnan(df["MS1"])
        tmt_channels = utils.filter_for_intensity_columns(df)
        df.loc[has_ms1, tmt_channels.columns] = tmt_channels.loc[has_ms1, :].multiply(
            df.loc[has_ms1, "MS1"] / tmt_channels.loc[has_ms1, :].sum(axis=1), axis="index"
        )

        return df
