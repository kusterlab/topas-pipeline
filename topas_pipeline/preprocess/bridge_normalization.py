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

from .. import sample_annotation
from . import phospho_grouping

# hacky way to get the package logger instead of just __main__ when running as a module
logger = logging.getLogger(__package__ + "." + Path(__file__).stem)


def apply_bridge_channel_normalization(
    results_folder: str,
    sample_annotation_file: str,
    min_occurrence: float = 2 / 3,  # Good Value for phospho
    min_samples_in_qc_lot: int = 8,
    overwrite: bool = False,
):
    results_folder = Path(results_folder)
    batch_corrected_file = results_folder / "preprocessed_pp2_agg_batchcorrected.csv"
    if batch_corrected_file.is_file():
        if not overwrite:
            logger.info(
                f"Reusing previously generated batch corrected intensities: {batch_corrected_file}"
            )
            return
        logger.info(f"Found existing results but overwrite flag was set.")

    phospho_df = phospho_grouping.read_cohort_intensities_df(
        f"{results_folder}/preprocessed_pp2_agg.csv"
    )

    sample_annotation_df = sample_annotation.load_sample_annotation(
        sample_annotation_file
    )
    sample_qc_lot_mapping_df = sample_annotation.get_sample_qc_lot_mapping_df(
        phospho_df.columns, sample_annotation_df
    )

    ref_cols = sample_qc_lot_mapping_df.loc[
        sample_qc_lot_mapping_df["is_reference"]
    ].index
    patient_columns = sample_qc_lot_mapping_df.loc[
        ~sample_qc_lot_mapping_df["is_reference"]
    ].index

    phospho_df_corrected, correction_factors = within_qc_lot_normalization(
        phospho_df, ref_cols, min_occurrence, sample_qc_lot_mapping_df
    )
    correction_factors.to_csv(results_folder / "peptide_correction_factors.csv")

    logger.info("Applying across QC lot normalization")

    # store original column order for output dataframe
    combined_column_order = phospho_df_corrected.columns

    # temporarily center patient data around 0 for each peptide
    avg_patient_intensity = phospho_df_corrected.loc[:, patient_columns].median(axis=1)
    phospho_df_corrected2 = (phospho_df_corrected.T - avg_patient_intensity).T

    # correct QC lots with at least 100 patient samples
    qc_lots_to_correct = (
        sample_qc_lot_mapping_df["QC Lot group"]
        .value_counts()[
            sample_qc_lot_mapping_df["QC Lot group"].value_counts()
            >= min_samples_in_qc_lot
        ]
        .index
    )
    for qc_lot_group in qc_lots_to_correct:
        qc_lot_samples = sample_qc_lot_mapping_df["QC Lot group"] == qc_lot_group
        qc_lot_patient_samples = sample_qc_lot_mapping_df.loc[
            qc_lot_samples & ~sample_qc_lot_mapping_df["is_reference"], :
        ].index

        avg_lot_patient_intensity = phospho_df_corrected2.loc[
            :, qc_lot_patient_samples
        ].median(axis=1)
        phospho_df_corrected2.loc[:, qc_lot_samples] = phospho_df_corrected2.loc[
            :, qc_lot_samples
        ].sub(avg_lot_patient_intensity, axis=0)

    # bring back to original intensity level
    phospho_df_corrected2 = (phospho_df_corrected2.T + avg_patient_intensity).T
    phospho_df_corrected2 = phospho_df_corrected2[combined_column_order]

    logger.info(f"Writing results to {batch_corrected_file}")
    phospho_df_corrected2.to_csv(
        batch_corrected_file,
        float_format="%.6f",
        index=True,
    )


def within_qc_lot_normalization(
    phospho_df: pd.DataFrame,
    ref_cols: list[str],
    min_occurrence: float,
    sample_qc_lot_mapping_df: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.Series]:
    """Filter for high-occurring p-peptides and apply row-wise normalization per batch and per QC lot.

    Keeps p-peptides that occur in at least min_occurrence of reference channels.

    Args:
        phospho_df (pd.DataFrame): _description_
        ref_cols (list[str]): _description_
        min_occurrence (float): _description_
        sample_qc_lot_mapping_df (pd.DataFrame): _description_

    Returns:
        tuple[pd.DataFrame, pd.Series]: _description_
    """
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

    # Calculates the ref correction for each mix type.
    # Aggregation happens as mean in the realspace. (Currently slow implementation !! TODO: pull out power and log)
    # Real space aggregation is more robust against low-intense outliers than mean in log space

    # Calculate the mean ref intensity per Batch and QC Lot
    mean_intensities = (
        vals[vals["is_reference"].astype("boolean").fillna(False)]
        .groupby(["QC Lot", "batch"], as_index=False)["intensity raw"]
        .agg(mean_realspace_f, engine="cython")
    )

    # Calculate the median ref intensity per QC lot
    mean_intensities["mix mean intensities"] = mean_intensities.groupby("QC Lot")[
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


def read_cohort_modified_sequence_groups(results_folder: Path) -> pd.DataFrame:
    df_patients = pd.read_csv(
        results_folder / "preprocessed_pp2_agg_batchcorrected.csv",
        usecols=["Gene names", "Modified sequence group"],
    )
    return df_patients[["Modified sequence group"]]


"""
python3 -m topas_pipeline.bridge_normalization -c config_patients.json
"""
if __name__ == "__main__":
    import argparse

    from .. import config

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-c", "--config", required=True, help="Absolute path to configuration file."
    )
    parser.add_argument(
        "-o",
        "--overwrite",
        action="store_true",
        help="Ignore existing results and recompute outputs.",
    )
    args = parser.parse_args(sys.argv[1:])

    configs = config.load(args.config)

    apply_bridge_channel_normalization(
        results_folder=configs.results_folder,
        sample_annotation_file=configs.sample_annotation,
        overwrite=args.overwrite,
    )
