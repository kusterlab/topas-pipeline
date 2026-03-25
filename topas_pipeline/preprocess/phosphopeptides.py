import time
import logging
from pathlib import Path

import pandas as pd

from . import bridge_normalization
from . import phospho_grouping
from .. import identification_metadata

# hacky way to get the package logger instead of just __main__ when running as a module
logger = logging.getLogger(__package__ + "." + Path(__file__).stem)


def update_with_batch_corrected_intensities(
    cohort_intensities_df: pd.DataFrame, cohort_batch_corrected_df: pd.DataFrame
) -> pd.DataFrame:
    return pd.concat(
        [
            cohort_intensities_df.loc[
                ~cohort_intensities_df.index.isin(cohort_batch_corrected_df.index)
            ],
            cohort_batch_corrected_df.loc[
                cohort_batch_corrected_df.index.isin(cohort_intensities_df.index)
            ].join(
                cohort_intensities_df.filter(
                    like=identification_metadata.METADATA_COLUMN_PREFIX
                ),
                how="left",
            ),
        ]
    )


def get_cohort_intensities_df(
    results_folder: str,
    sample_annotation_file: str = None,
    keep_identification_metadata_columns: bool = True,
):
    cohort_intensities_df = phospho_grouping.read_cohort_intensities_df(
        f"{results_folder}/preprocessed_pp2_agg.csv",
        sample_annotation_file,
        keep_identification_metadata_columns=keep_identification_metadata_columns,
    )

    cohort_batch_corrected_df = phospho_grouping.read_cohort_intensities_df(
        f"{results_folder}/preprocessed_pp2_agg_batchcorrected.csv",
        sample_annotation_file,
        keep_identification_metadata_columns=False,
    )

    cohort_intensities_df = update_with_batch_corrected_intensities(
        cohort_intensities_df, cohort_batch_corrected_df
    )

    return cohort_intensities_df


def group_phosphopeptides_and_normalize(
    results_folder: str,
    sample_annotation_file: str,
    run_lfq: bool,
):
    """
    extra phospho processing (grouping, normalization, expression correction)
    """
    # ~15 minutes for 2000 samples
    start_time = time.time()
    phospho_grouping.aggregate_modified_sequences(results_folder=results_folder)
    logger.info("--- %.1f seconds --- phospho grouping" % (time.time() - start_time))

    if run_lfq:
        start_time = time.time()
        filter_for_occurrence(
            results_folder=results_folder,
        )
        logger.info(
            "--- %.1f seconds --- filter for occurrence" % (time.time() - start_time)
        )
    else:
        start_time = time.time()
        bridge_normalization.apply_bridge_channel_normalization(
            results_folder=results_folder,
            sample_annotation_file=sample_annotation_file,
        )
        logger.info(
            "--- %.1f seconds --- bridge normalization" % (time.time() - start_time)
        )

    start_time = time.time()
    df = get_cohort_intensities_df(results_folder, sample_annotation_file)
    logger.info(
        "--- %.1f seconds --- combining non-normalized and bridge normalized p-peptide groups"
        % (time.time() - start_time)
    )

    return df


def filter_for_occurrence(
    results_folder: str,
    min_occurrence: float = 2 / 3,  # Good Value for phospho
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
    high_count_mask = phospho_df.notnull().mean(axis=1) > min_occurrence
    phospho_df = phospho_df[high_count_mask]

    logger.info(f"Writing results to {batch_corrected_file}")
    phospho_df.to_csv(
        batch_corrected_file,
        float_format="%.6f",
        index=True,
    )
