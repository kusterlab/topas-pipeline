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
    sample_annotation_file: str,
    keep_identification_metadata_columns: bool = True,
):
    cohort_intensities_df = phospho_grouping.read_cohort_intensities_df(
        f"{results_folder}/preprocessed_pp2_agg.csv",
        sample_annotation_file,
        keep_identification_metadata_columns=keep_identification_metadata_columns,
    )

    cohort_batch_corrected_df = phospho_grouping.read_cohort_intensities_df(
        f"{results_folder}/preprocessed_pp2_agg_batchcorrected.csv",
        keep_identification_metadata_columns=False,
    )

    cohort_intensities_df = update_with_batch_corrected_intensities(
        cohort_intensities_df, cohort_batch_corrected_df
    )

    return cohort_intensities_df


def group_phosphopeptides_and_normalize(
    results_folder: str, sample_annotation_file: str
):
    """
    extra phospho processing (grouping, normalization, expression correction)
    """
    # ~15 minutes for 2000 samples
    start_time = time.time()
    phospho_grouping.aggregate_modified_sequences(results_folder=results_folder)
    logger.info("--- %.1f seconds --- phospho grouping" % (time.time() - start_time))

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
