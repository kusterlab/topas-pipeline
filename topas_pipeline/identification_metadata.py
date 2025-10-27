import re
import logging

import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression

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
    detected_in_batch_df = detected_in_batch_df.fillna(False)
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


def mark_quant_out_of_range(df, ratio_threshold: float = 0.01):
    """
    Annotate proteins with a log ratio below the specified threshold.
    
    Parameters:
    df (DataFrame): The input DataFrame containing evidence data.
    logratio (int): The log ratio threshold for filtering.
    
    Returns:
    DataFrame: Filtered DataFrame with high DR evidence.
    """

    logger.info("Mark proteins with high dynamic range in the batch (lower than threshold)")

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
        tmt_cols_df = 10 ** tmt_cols_df

        metadata_cols = get_metadata_columns(df.filter(regex=rf"{batch}$"))
        
        row_max = tmt_cols_df.max(axis=1)
        ratio_df = tmt_cols_df.div(row_max, axis=0).astype(float)
        mask = ratio_df < ratio_threshold

        df[metadata_cols.columns] += np.where(
            mask, "quan_OOR;", ""
        )

    return df


def mark_peptide_id_out_of_range(df, threshold: int = 2):

    logger.info("Mark proteins with peptide IDs out of range")
    
    metadata_cols_df = get_metadata_columns(df)
    n_meta_cols = len(metadata_cols_df.columns)

    pep_id_pattern = r'num_peptides=(\d+)'
    num_peptides = metadata_cols_df.copy()
    num_peptides = num_peptides.astype(str).applymap(
        lambda x: re.search(pep_id_pattern, x).group(1) if re.search(pep_id_pattern, x) else None
    )
    df['is_outlier_peptide_id'] = None

    # this would be the function  - give df, num_peptides, threshold
    for i, row in df.iterrows(): 
        y = row.loc[(row.index.str.contains(r'^Reporter'))].astype(float).values
        y = np.where(y == 0, np.nan, y)
        X = num_peptides.loc[i, :].values.reshape(-1, 1)

        # create a mask for non-nan values
        mask = ~np.isnan(y)

        # filter out NaN in both
        X_sub = X[mask].astype(int)
        y_sub = y[mask]

        # Skip if insufficient data or no variation
        if len(y_sub) < 5 or np.unique(X_sub).size == 1:
            df.at[i, 'is_outlier_peptide_id'] = [False] * n_meta_cols
            continue
        
        # fit regression
        linear = LinearRegression()
        linear.fit(X_sub, y_sub)

        # predict and calculate residuals
        predicted = linear.predict(X_sub)
        residuals = y_sub - predicted


        # standardize residuals
        residual_std = np.nanstd(residuals)
        if residual_std == 0 or np.isnan(residual_std):
            df.at[i, 'is_outlier_peptide_id'] = [False] * n_meta_cols
            continue
        standardized_residuals = residuals / residual_std

        # weighted threshold
        weighted_thresholds = threshold * (1 + (np.power(1.1, X_sub - 1) - 1))

        # flag outliers
        is_outlier = standardized_residuals > weighted_thresholds.ravel()

        full_outliers = np.full(n_meta_cols, False)
        full_outliers[mask] = is_outlier

        df.at[i, 'is_outlier_peptide_id'] = full_outliers.tolist()


    # instead make a unit test for this
    # assert all(df['is_outlier_peptide_id'].apply(len) == n_meta_cols), \
    # "Each boolean list must match number of metadata columns"

    # Expand the boolean lists into a DataFrame aligned by index
    expanded = pd.DataFrame(df['is_outlier_peptide_id'].tolist(),
                            index=df.index,
                            columns=metadata_cols_df.columns).fillna(False).astype(bool)

    # Append "quan_OOR;" only where True
    df.loc[:, metadata_cols_df.columns] = (
        df.loc[:, metadata_cols_df.columns]
        .mask(expanded, df[metadata_cols_df.columns].astype(str) + "pepcount_OOR;")
    )

    return df
