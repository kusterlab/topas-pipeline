from abc import ABC, abstractmethod
import pandas as pd

from topas_pipeline.preprocess.preprocess_tools import (
    merge_ms1_columns,
    median_centering,
)
from topas_pipeline.utils import filter_for_intensity_columns


class BaseQuantDataLoader(ABC):

    @abstractmethod
    def load_and_normalize(self) -> pd.DataFrame:
        raise NotImplementedError("Subclasses must implement this method")

    @abstractmethod
    def load(self) -> pd.DataFrame:
        raise NotImplementedError("Subclasses must implement this method")

    @abstractmethod
    def normalize(
        self,
        df: pd.DataFrame,
        ref_channels_df: pd.DataFrame,
        normalize_to_reference: bool,
    ) -> pd.DataFrame:
        raise NotImplementedError("Subclasses must implement this method")

    def median_centering_within_batch(self, df_raw: pd.DataFrame) -> pd.DataFrame:
        """Normalizes samples by median centering within each batch.
        For each batch, calculates the median of each TMT channel and computes a correction factor
        to center the median of each channel to the overall median across all channels.
        Applies the correction factor to the TMT channel intensities within each batch.
        """
        dfs = []

        # Get column names that contain TMT channel intensities
        tmt_channels = filter_for_intensity_columns(df_raw).columns

        correction_factors_all = {}
        for batch_name, df in df_raw.groupby("Batch"):
            df.loc[:, tmt_channels], correction_factors = median_centering(
                df.loc[:, tmt_channels]
            )
            correction_factors_all[batch_name] = correction_factors
            dfs.append(df)

        df_normalized = pd.concat(dfs, axis=0)

        correction_factors_df = pd.DataFrame.from_dict(correction_factors_all)
        correction_factors_df = correction_factors_df.reset_index()
        correction_factors_df = pd.melt(
            correction_factors_df,
            id_vars="index",
            var_name="Variable",
            value_name="Value",
        )

        return df_normalized, correction_factors_df

    def median_centering_ms1(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Normalizes samples by multiplying each batch MS1 with its own correction factor.
        Only uses peptides detected in >70% of the batches for computing median to prevent
        low-abundant peptides from dragging down the median in runs with deeper coverage
        :param all_batches: list of evidence dataframes
        :return:
        """

        merged_ms1_df = merge_ms1_columns(df)

        medians = (
            merged_ms1_df[
                merged_ms1_df.count(axis=1) > 0.7 * len(merged_ms1_df.columns)
            ]
            .median(axis=0)
            .to_dict()
        )
        mean_median = pd.Series(medians.values()).mean()

        dfs, correction_factors = [], {}
        for _, df in df.groupby("Batch"):
            batch_name = df["Batch"].iloc[0]
            correction_factor = mean_median / medians[batch_name]

            df["MS1"] = df["Intensity"] * correction_factor
            dfs.append(df)
            correction_factors[batch_name] = correction_factor

        df_normalized = pd.concat(dfs, axis=0)

        correction_factors_df = pd.DataFrame(correction_factors, index=[0]).T
        return df_normalized, correction_factors_df

    def scale_ms1_to_reference(
        self, df: pd.DataFrame, ref_channel_df: pd.DataFrame
    ) -> pd.DataFrame:
        return df

    def impute_ms1_intensity(
        self, df: pd.DataFrame, ref_channel_df: pd.DataFrame
    ) -> pd.DataFrame:
        return df

    def redistribute_ms1_intensity(self, df: pd.DataFrame) -> pd.DataFrame:
        return df
