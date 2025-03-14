import re
import logging

import pandas as pd
import numpy as np

from . import utils

logger = logging.getLogger(__name__)
METADATA_COLUMN_PREFIX = "Identification metadata"


def create_metadata_columns(df: pd.DataFrame) -> pd.DataFrame:
    """create columns to store metadata about the identifications, e.g. imputed, detected in batch, single peptide id
    annotations should be separated by semi-colons, e.g. imputed;single peptide id;detected in batch;
    """
    metadata_cols = as_metadata_columns(utils.get_tmt_channels(df).columns.str)
    df[metadata_cols] = ""
    return df


def as_metadata_columns(x: str):
    return x.replace("Reporter intensity corrected", METADATA_COLUMN_PREFIX)


def remove_metadata_column_prefix(df: pd.DataFrame) -> pd.DataFrame:
    return df.rename(columns=lambda x: x.replace(f"{METADATA_COLUMN_PREFIX} ", ""))


def get_metadata_columns(df: pd.DataFrame) -> pd.DataFrame:
    metadata_df = df.filter(regex=rf"^{METADATA_COLUMN_PREFIX}")
    if len(metadata_df.columns) == 0:
        raise ValueError(
            f"Could not find '{METADATA_COLUMN_PREFIX}' columns in the dataframe"
        )
    return metadata_df


def mark_detected_in_batch(df: pd.DataFrame) -> pd.DataFrame:
    logger.info("Mark proteins detected in the batch but not in the sample")
    batches = {
        re.sub(r"Reporter intensity corrected \d{1,2} ", "", c)
        for c in df.columns
        if c.startswith("Reporter intensity corrected")
    }

    for batch in batches:
        tmt_cols_df = df.filter(
            regex=rf"^Reporter intensity corrected \d{{1,2}} {batch}$", axis="columns"
        )
        tmt_cols_df = tmt_cols_df.replace(0, np.nan)

        num_detected_in_batch = tmt_cols_df.count(axis="columns")
        imputation_cols = as_metadata_columns(tmt_cols_df.columns.str)
        df[imputation_cols] += np.where(
            tmt_cols_df.isna().multiply(num_detected_in_batch > 0, axis="index"),
            "detected in batch;",
            "",
        )

    return df


def get_detected_in_batch(df: pd.DataFrame) -> pd.DataFrame:
    metadata_df = get_metadata_columns(df)
    return metadata_df.apply(lambda x: x.str.contains("detected in batch;"), axis=1)


def replace_detected_in_batch(
    df: pd.DataFrame, annot_df: pd.DataFrame, value: float
) -> pd.DataFrame:
    """Replace nan values for peptides/proteins detected in a batch by a value"""
    detected_in_batch_df = get_detected_in_batch(annot_df)
    detected_in_batch_df = remove_metadata_column_prefix(detected_in_batch_df)

    patient_columns_prefix = detected_in_batch_df.add_prefix("pat_")
    detected_in_batch_df = detected_in_batch_df.rename(
        dict(zip(detected_in_batch_df.columns, patient_columns_prefix)), axis=1
    )

    detected_in_batch_df = detected_in_batch_df[df.columns].astype('object')
    detected_in_batch_df[:] = np.where(detected_in_batch_df, value, np.nan)
    df.update(
        detected_in_batch_df, overwrite=False
    )  # for some reason update can only be done in place
    return df


def mark_num_peptides(df: pd.DataFrame) -> pd.DataFrame:
    logger.info("Mark number of peptide identifications per sample")
    peptide_counts_df = df.filter(
        regex=rf"^Unique peptides Reporter intensity corrected \d{{1,2}}",
        axis="columns",
    )
    imputation_cols = as_metadata_columns(
        peptide_counts_df.columns.str.replace("Unique peptides ", "").str
    )

    def num_peptide_string(peptide_count_cols):
        return [f"num_peptides={int(c)};" if c > 0 else "" for c in peptide_count_cols]

    df[imputation_cols] += peptide_counts_df.apply(
        num_peptide_string, axis=1, result_type="expand"
    ).values
    return df


def get_num_peptides(df: pd.DataFrame) -> pd.DataFrame:
    metadata_df = get_metadata_columns(df)
    return metadata_df.apply(
        lambda x: x.str.extract(r"num_peptides=(\d+);")[0], axis=1
    ).astype(float)


def filter_by_min_peptides(
    df: pd.DataFrame, annot_df: pd.DataFrame, min_peptides: int
) -> pd.DataFrame:
    """Replace nan values for peptides/proteins detected in a batch by a value"""
    detected_in_batch_df = get_num_peptides(annot_df)
    detected_in_batch_df = remove_metadata_column_prefix(detected_in_batch_df)
    detected_in_batch_df = detected_in_batch_df[df.columns]
    detected_in_batch_df = detected_in_batch_df.applymap(
        lambda x: -1 if x < min_peptides else np.nan
    )
    df.update(
        detected_in_batch_df, overwrite=True
    )  # for some reason update can only be done in place
    df = df.replace(-1, np.nan)
    return df
