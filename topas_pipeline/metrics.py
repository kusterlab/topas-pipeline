import os.path
import os
import sys
import warnings
import logging
from pathlib import Path
from typing import List, Dict, Tuple, Union

import scipy.stats
import numpy as np
import pandas as pd

from . import clinical_tools
from . import clinical_annotation
from . import utils

logger = logging.getLogger(__name__)

MEASURE_NAMES = [
    ("rank", "rank"),
    ("fc", "fold_change"),
    ("z", "z-score"),
    ("p", "p-value"),
]


def compute_metrics(
    results_folder: Union[str, Path], debug: bool, data_types: List[str]
):
    """Computes rank, fold change, z-score and p-value
    Requires that annot_{data_type}.csv is generated by clinical_annotation.py
    """
    logger.info("Running compute metrics module")

    if debug:
        data_types += [data_type + "_with_ref" for data_type in data_types]

    for data_type in data_types:
        if check_measures_computed(results_folder, data_type, MEASURE_NAMES):
            logger.info(f"Skipping compute_metrics {data_type} - already computed")
            continue

        logger.info(f"Reading in clinically annotated {data_type} data")
        annot_df = clinical_annotation.read_annotated_expression_file(
            results_folder, data_type
        )

        measures = get_metrics(annot_df.filter(regex=r"(^pat_)|(^ref_)"))
        # TODO: check if these annotation columns are still used anywhere, they 
        # have been absent since renaming 'BASKET' to 'TOPAS_SCORE' in the TOPAS
        # annotation file around January 2025
        # measures["z-score"] = add_topas_annotations(measures["z-score"], annot_df)

        save_measures(
            results_folder,
            MEASURE_NAMES,
            measures,
            data_type,
        )


def add_topas_annotations(
    measure_df: pd.DataFrame, annot_df: pd.DataFrame
) -> pd.DataFrame:
    topas_score_col = clinical_tools.TOPAS_SCORE_COLUMN
    topas_subscore_col = clinical_tools.TOPAS_SUBSCORE_COLUMN
    measure_df = measure_df.join(
        annot_df.loc[
            :,
            [
                topas_score_col,
                f"{topas_score_col}_weights",
                topas_subscore_col,
                f"{topas_subscore_col}_weights",
            ],
        ]
    )
    return measure_df


class Metrics:
    """
    Describe
    """

    def __init__(self, df):
        self.df = df
        self.metrics_df = {}
        self.calc()

    @staticmethod
    def get_rank(df: pd.DataFrame) -> pd.DataFrame:
        """ """
        logger.debug("Calculating ranks")
        df_rank = df.copy()
        df_rank = df_rank.rank(
            ascending=False, method="average", na_option="keep", axis=1
        )
        # add new column that for each protein/peptide tells the max rank
        df_rank["max"] = df_rank.notnull().sum(axis=1)
        return df_rank.add_prefix("rank_")

    @staticmethod
    def get_fold_change(df: pd.DataFrame) -> pd.DataFrame:
        """
        Leave-one-out approach (loo)
        """
        logger.debug("Calculating fold changes")
        df_fold_change = np.power(10, df)
        df_fold_change /= _loo_median(df_fold_change.values)
        df_fold_change = df_fold_change.add_prefix("fc_")
        return df_fold_change

    @staticmethod
    def get_zscore(df: pd.DataFrame) -> pd.DataFrame:
        """
        Leave-one-out approach (loo)
        """
        logger.debug("Calculating z-scores")
        df_z_score = df.copy()
        df_z_score -= _loo_median(df.values)
        df_z_score /= _loo_std(df.values)
        df_z_score = df_z_score.add_prefix("zscore_")
        return df_z_score

    @staticmethod
    def get_pvalues(z_scores: pd.DataFrame) -> pd.DataFrame:
        """ """
        logger.debug("Calculating p-values")
        df_p_values = scipy.stats.norm.sf(abs(z_scores)) * 2
        df_p_values = pd.DataFrame(df_p_values)
        df_p_values.columns = z_scores.columns
        df_p_values.index = z_scores.index
        return df_p_values.add_prefix("pvalue_")

    def calc(self):
        # add metrics
        self.metrics_df["rank"] = Metrics.get_rank(self.df)
        self.metrics_df["fold_change"] = Metrics.get_fold_change(self.df)
        self.metrics_df["z-score"] = Metrics.get_zscore(self.df)
        self.metrics_df["p-value"] = Metrics.get_pvalues(self.metrics_df["z-score"])
        self.metrics_df["occurrence"] = self.df.count(axis=1)
        self.metrics_df["occurrence"].name = "occurrence"
        self.metrics_df["occurrence"] = self.metrics_df["occurrence"].to_frame()


def _loo_median(input_matrix: np.array) -> np.array:
    """
    Computes the leave-one-out median for each element row-wise.

    30 seconds for 100000 proteins and 2000 samples.
    """
    # 1. get the rank of each value within the row
    argsort_indices = np.argsort(input_matrix, axis=1)
    ranks = np.argsort(argsort_indices, axis=1)

    valid_vals = ranks.shape[1] - np.isnan(input_matrix).sum(axis=1)[:, np.newaxis]

    # 2. get the LOO median indices
    # lower_median and upper_median indices are the same if the LOO array has an odd number of elements
    lower_median = np.floor(
        valid_vals / 2  # index of median
        - (
            ranks >= valid_vals / 2
        )  # if the value is in the upper half, go one element down
        - (
            (ranks * 2 + 1) == valid_vals
        )  # if the value is exactly the median, go another element down
    ).astype(int)
    upper_median = np.ceil(valid_vals / 2 - (ranks >= valid_vals / 2)).astype(int)

    sorted_values = input_matrix[
        np.arange(input_matrix.shape[0])[:, np.newaxis], argsort_indices
    ]

    # reuse variable to (hopefully) save memory
    loo_median = sorted_values[
        np.arange(input_matrix.shape[0])[:, np.newaxis], lower_median
    ]
    loo_median += sorted_values[
        np.arange(input_matrix.shape[0])[:, np.newaxis], upper_median
    ]
    loo_median /= 2
    loo_median[np.isnan(input_matrix)] = np.nan

    return loo_median


def _loo_std(input_matrix: np.array) -> np.array:
    """LOO standard deviation using Welford's online algorithm reversed

    Args:
        input_matrix (np.array): _description_

    Returns:
        _type_: _description_
    """
    n = input_matrix.shape[1] - np.isnan(input_matrix).sum(axis=1)[:, np.newaxis]

    regular_sum = np.nansum(input_matrix, axis=1)[:, np.newaxis]
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        regular_var = np.nanvar(input_matrix, axis=1)[
            :, np.newaxis
        ]  # uses ddof = 0 because that is the numpy default
    m2n = regular_var * n

    loo_sum = regular_sum - input_matrix
    loo_mean = loo_sum / np.clip(n - 1, 1, None)
    m2n_minus_1 = m2n - (input_matrix - loo_mean) * (
        input_matrix - regular_sum / np.clip(n, 1, None)
    )

    return np.sqrt(
        m2n_minus_1 / np.where(n - 2 >= 1, n - 2, np.nan)
    )  # uses ddof = 1 because that is the pandas default


def get_metrics(df: pd.DataFrame) -> Dict[str, pd.DataFrame]:
    """
    describe
    """
    logger.info("Calculating metrics")
    df = utils.keep_only_sample_columns(df)

    # Get metrics
    m = Metrics(df)
    return m.metrics_df


def get_data_type_long(data_type: str):
    # TODO: raise error for unknown data type instead of defaulting to full_proteome
    if data_type.startswith("pp"):
        return "phospho"
    return "full_proteome"


def check_measures_computed(
    results_folder: Union[str, Path],
    data_type: str,
    measure_names: List[Tuple[str]] = MEASURE_NAMES,
) -> Dict[str, pd.DataFrame]:
    """
    explain what modalities and measures have to be

    """
    data_type_long = get_data_type_long(data_type)
    for m, _ in measure_names:
        filename = os.path.join(results_folder, f"{data_type_long}_measures_{m}.tsv")
        if not os.path.exists(filename):
            return False

    return True


def read_measures(
    results_folder: Union[str, Path],
    data_type: str,
    measure_names: List[Tuple[str]] = MEASURE_NAMES,
    with_reference_channels: bool = False,
) -> Dict[str, pd.DataFrame]:
    """
    explain what modalities and measures have to be

    """
    data_type_long = get_data_type_long(data_type)
    index_col = utils.get_index_cols(data_type)

    measures = dict()
    for m, measure in measure_names:
        filename = os.path.join(results_folder, f"{data_type_long}_measures_{m}.tsv")
        if not os.path.exists(filename):
            return dict()

        logger.info(f"Reading in {filename}")
        measures[measure] = pd.read_csv(filename, sep="\t", index_col=index_col)
    return measures


def get_topas_annotation_columns() -> List[str]:
    topas_score_col = list(clinical_tools.TOPAS_SCORE_COLUMNS.keys())[0]
    topas_subscore_col = list(clinical_tools.TOPAS_SUBSCORE_COLUMNS.keys())[0]
    return [
        topas_score_col,
        f"{topas_score_col}_weights",
        topas_subscore_col,
        f"{topas_subscore_col}_weights",
    ]


def save_measures(
    results_folder: Union[str, Path],
    measure_names: List[Tuple[str]],
    measures: Dict[str, pd.DataFrame],
    data_type: str,
):
    # save quantification measures
    data_type_long = get_data_type_long(data_type)

    filename_suffix = ""
    if data_type.endswith("_with_ref"):
        filename_suffix = "_ref"

    for m, measure in measure_names:
        filename = os.path.join(
            results_folder, f"{data_type_long}_measures_{m}{filename_suffix}.tsv"
        )
        measures[measure].to_csv(filename, sep="\t", float_format="%.6g")


if __name__ == "__main__":
    import argparse

    from . import config

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-c", "--config", required=True, help="Absolute path to configuration file."
    )
    args = parser.parse_args(sys.argv[1:])

    configs = config.load(args.config)

    compute_metrics(
        configs.results_folder,
        configs.preprocessing.debug,
        data_types=configs.data_types,
    )
