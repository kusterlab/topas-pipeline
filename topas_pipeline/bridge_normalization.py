import sys
import logging
from pathlib import Path

import pandas as pd
import numpy as np

from tqdm import tqdm

tqdm.pandas()

import warnings

warnings.simplefilter(action="ignore", category=FutureWarning)
warnings.filterwarnings(action="ignore", category=RuntimeWarning)

from . import sample_annotation
from . import preprocess_tools as prep
from . import phospho_grouping

# hacky way to get the package logger instead of just __main__ when running as a module
logger = logging.getLogger(__package__ + "." + Path(__file__).stem)


def apply_bridge_channel_normalization(
    results_folder: str,
    sample_annotation_file: str,
    min_occurrence: float = 2 / 3,  # Good Value for phospho
):
    results_folder = Path(results_folder)
    sample_annotation_df = sample_annotation.load_sample_annotation(
        sample_annotation_file
    )

    phospho_df = phospho_grouping.read_cohort_intensities_df(results_folder)

    sample_qc_lot_mapping_df = get_sample_qc_lot_mapping_df(
        phospho_df.columns, sample_annotation_df
    )

    ref_cols = sample_qc_lot_mapping_df.loc[
        sample_qc_lot_mapping_df["QC Lot"].notna()
    ].index
    pat_cols = sample_qc_lot_mapping_df.loc[
        sample_qc_lot_mapping_df["QC Lot"].isna()
    ].index

    phospho_df_corrected, correction_factors = within_qc_lot_normalization(
        phospho_df, ref_cols, min_occurrence, sample_qc_lot_mapping_df
    )
    correction_factors.to_csv(results_folder / "peptide_correction_factors.csv")

    logger.info("Applying across QC lot normalization")

    combined_column_order = phospho_df_corrected.columns

    phospho_df_ref_corrected = phospho_df_corrected.loc[:, ref_cols]
    phospho_df_corrected = phospho_df_corrected.loc[:, pat_cols]

    # temporarily center patient data around 0 for each peptide
    avg_intensity = phospho_df_corrected.median(axis=1)
    phospho_df_corrected2 = (phospho_df_corrected.T - avg_intensity).T
    phospho_df_ref_corrected2 = (phospho_df_ref_corrected.T - avg_intensity).T

    # correct QC lots with at least 100 patient samples
    qc_lots_to_correct = (
        sample_qc_lot_mapping_df["QC Lot group"]
        .value_counts()[sample_qc_lot_mapping_df["QC Lot group"].value_counts() >= 100]
        .index
    )
    for qc_lot_group in qc_lots_to_correct:
        qc_lot_samples = sample_qc_lot_mapping_df["QC Lot group"] == qc_lot_group
        qc_lot_patient_samples = sample_qc_lot_mapping_df.loc[
            qc_lot_samples & sample_qc_lot_mapping_df["QC Lot"].isna()
        ].index
        qc_lot_ref_samples = sample_qc_lot_mapping_df.loc[
            qc_lot_samples & sample_qc_lot_mapping_df["QC Lot"].notna()
        ].index

        avg_patient_intensity = phospho_df_corrected2.loc[
            :, qc_lot_patient_samples
        ].median(axis=1)
        phospho_df_corrected2.loc[:, qc_lot_patient_samples] = (
            phospho_df_corrected2.loc[:, qc_lot_patient_samples].sub(
                avg_patient_intensity, axis=0
            )
        )
        phospho_df_ref_corrected2.loc[:, qc_lot_ref_samples] = (
            phospho_df_ref_corrected2.loc[:, qc_lot_ref_samples].sub(
                avg_patient_intensity, axis=0
            )
        )

    # bring back to original intensity level
    phospho_df_corrected2 = (phospho_df_corrected2.T + avg_intensity).T
    phospho_df_ref_corrected2 = (phospho_df_ref_corrected2.T + avg_intensity).T

    phospho_df_corrected2 = pd.concat(
        [phospho_df_corrected2, phospho_df_ref_corrected2], axis=1
    )
    phospho_df_corrected2 = phospho_df_corrected2[combined_column_order]

    # map column names to patient names
    channel_to_sample_id_dict = sample_annotation.get_channel_to_sample_id_dict(
        sample_annotation_df,
        remove_qc_failed=True,
        remove_replicates=False,
    )

    index_cols = ["Modified sequence group", "Gene names"]
    phospho_df_corrected2 = prep.rename_columns_with_sample_ids(
        phospho_df_corrected2.reset_index(),
        channel_to_sample_id_dict,
        index_cols=index_cols,
    )

    batch_corrected_file = results_folder / "preprocessed_pp2_agg_batchcorrected.csv"
    logger.info(f"Writing results to {batch_corrected_file}")
    phospho_df_corrected2.to_csv(
        batch_corrected_file,
        float_format="%.6f",
        index=False,
    )


def get_sample_qc_lot_mapping_df(
    sample_columns: pd.Index, sample_annotation_df: pd.DataFrame
) -> pd.DataFrame:
    sample_mapping_df = sample_columns.to_frame().reset_index()
    sample_mapping_df["batch"] = sample_mapping_df["index"].str.split("_Batch").str[-1]
    sample_mapping_df["channel"] = (
        sample_mapping_df["index"].str.split(" ").str[-2].astype(int)
    )
    sample_mapping_df = sample_mapping_df.merge(
        sample_annotation_df[["Batch Name", "TMT Channel", "QC Lot"]],
        left_on=["batch", "channel"],
        right_on=["Batch Name", "TMT Channel"],
    )
    sample_mapping_df["QC Lot group"] = sample_mapping_df.groupby("batch")[
        "QC Lot"
    ].transform("mean")
    sample_mapping_df = sample_mapping_df.drop(columns=[0, "Batch Name", "TMT Channel"])
    sample_mapping_df = sample_mapping_df.set_index("index")
    return sample_mapping_df


def within_qc_lot_normalization(
    phospho_df: pd.DataFrame,
    ref_cols: list[str],
    min_occurrence: float,
    sample_qc_lot_mapping_df: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.Series]:
    logger.info("Applying within QC lot normalization")
    high_count_mask = phospho_df[ref_cols].notnull().mean(axis=1) > min_occurrence
    phospho_df_corrected = phospho_df[high_count_mask].progress_apply(
        row_wise_normalize, sample_qc_lot_mapping_df=sample_qc_lot_mapping_df, axis=1
    )
    correction_factors = phospho_df_corrected - phospho_df[high_count_mask]
    return phospho_df_corrected, correction_factors


def mean_f(x):
    return np.nanmean(x)


def mean_realspace_f(x):
    return np.log10(np.nanmean(np.power(10, x)))


def row_wise_normalize(
    row: pd.Series, sample_qc_lot_mapping_df: pd.DataFrame
) -> pd.Series:
    """
    return values : "intensity raw", "ref correction", "intensity corrected"
    """
    # Grab a single peptide
    vals = row.rename("intensity raw")
    vals = pd.concat([vals, sample_qc_lot_mapping_df], axis=1)
    vals["sample type"] = vals["QC Lot"]

    # Calculates the ref correction for each mix type.
    # Aggregation happens as mean in the realspace. (Currently slow implementation !! TODO: pull out power and log)
    # Real space aggregation is more robust against low-intense outliers than mean in log space

    # Calculate the mean ref intensity per Batch and QC Lot
    mean_intensities = (
        vals[vals["channel"].isin({10, 11})]
        .groupby(["sample type", "batch"], as_index=False)["intensity raw"]
        .agg(mean_realspace_f, engine="cython")
    )

    # Calculate the mean ref intensity per QC lot
    mean_intensities["mix mean intensities"] = mean_intensities.groupby("sample type")[
        "intensity raw"
    ].transform("median")

    # Calculate correction factor relative to mean of QC lot
    mean_intensities["ref correction"] = (
        mean_intensities["mix mean intensities"] - mean_intensities["intensity raw"]
    )

    # Aggregate multiple mixes into a single batch correction value (necessary if multiple QC lots are present in a single batch)
    batch_correction = mean_intensities.groupby("batch")["ref correction"].mean()
    vals = vals.merge(batch_correction, left_on="batch", right_index=True)

    # Correct intensity values
    vals["intensity corrected"] = vals["intensity raw"] + vals["ref correction"]
    return vals["intensity corrected"]


def read_cohort_batch_corrected_df(results_folder: str, skiprows: pd.Series = None):
    phospho_batch_corrected = pd.read_csv(
        results_folder / "preprocessed_pp2_agg_batchcorrected.csv",
        index_col=[0, 1],
        skiprows=skiprows,
    )
    return phospho_batch_corrected


"""
python3 -m topas_pipeline.bridge_normalization -c config_patients.json
"""
if __name__ == "__main__":
    import argparse

    from . import config

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-c", "--config", required=True, help="Absolute path to configuration file."
    )
    args = parser.parse_args(sys.argv[1:])

    configs = config.load(args.config)

    apply_bridge_channel_normalization(
        results_folder=configs.results_folder, sample_annotation_file=configs.sample_annotation
    )
