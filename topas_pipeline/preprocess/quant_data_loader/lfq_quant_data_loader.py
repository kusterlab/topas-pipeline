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


class LFQQuantDataLoader(BaseQuantDataLoader):
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

        # Get Result File loader based on quant file format used (Maxquant, Ionquant, Diann)
        result_file_loader_class = ResultFileLoaderFactory.get_loader(
            self.quant_file_format
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
        all_batches = self.result_file_loader.load()
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
