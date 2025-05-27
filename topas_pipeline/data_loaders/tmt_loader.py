from typing import Union, List
from pathlib import Path
import re
import logging

import numpy as np
import pandas as pd

from .data_loader import DataLoader, extract_batch_name
from ..utils import get_tmt_channels, get_ref_channels

logger = logging.getLogger(__name__)


class TMTLoader(DataLoader):
    def __init__(self, evidence_files):
        self.evidence_files = evidence_files

    def load_data(self, use_cols: List[str]):
        all_batches = [_parse_evidence_file(f, use_cols) for f in self.evidence_files]
        return all_batches

    def scale_ms1_to_reference(
        self, df: pd.DataFrame, ref_channel_df: pd.DataFrame
    ) -> pd.DataFrame:
        return _scale_ms1_to_reference(df, ref_channel_df)

    def impute_ms1_intensity(
        self, df: pd.DataFrame, ref_channel_df: pd.DataFrame
    ) -> pd.DataFrame:
        return _impute_ms1_intensity(df, ref_channel_df)

    def redistribute_ms1_intensity(self, df: pd.DataFrame) -> pd.DataFrame:
        return _redistribute_ms1_intensity(df)


def _parse_evidence_file(
    evidence_file: Union[str, Path], use_cols: List[str]
) -> pd.DataFrame:
    logger.info(f"Parsing evidence file: {evidence_file}")
    df = pd.read_csv(evidence_file, sep="\t", usecols=use_cols)

    df = df.astype(
        {
            "Reverse": "category",
            "Experiment": "category",
            "Modifications": "category",
            "Potential contaminant": "category",
        }
    )

    # Change phosphosite notation in modified sequence
    df["Modified sequence"] = df["Modified sequence"].str.replace(
        re.compile(r"([STY])\(Phospho \(STY\)\)"),
        lambda pat: f"p{pat.group(1)}",
        regex=True,
    )

    # Remove rows with no measured reporter ions
    df = df.replace(0, np.nan)
    df = df.loc[
        ~df.filter(regex="^Reporter intensity corrected").isnull().all(axis=1), :
    ]

    batch_name = extract_batch_name(evidence_file)
    df["Batch"] = batch_name
    return df


def _scale_ms1_to_reference(
    df: pd.DataFrame, ref_channel_df: pd.DataFrame
) -> pd.DataFrame:
    """Scale MS1 intensities such that reference channel MS1 contribution is constant across batches.

    TODO: combine this with _impute_ms1_intensity.

    For each modified peptide in each batch check how many batches have MS1 and reference channels.
    For those, compute the mean intensity of the reference channels after MS1 redistribution.
    :param df: dataframe of batch intensities (11 reporter channels, intensities are not log transformed and missing values as NaNs
                                               MS1 column, intensities are not log transformed and missing values as zeroes)
    """
    logger.info("Scale MS1 intensities using reference channels")
    df["MS1"] = df["MS1"].replace(0, np.nan)

    tmt_channels = get_tmt_channels(df)
    ref_channels = get_ref_channels(df, ref_channel_df)

    df["MS1 corrected reference intensity"] = (
        ref_channels.mean(axis=1) / tmt_channels.sum(axis=1) * df["MS1"]
    )

    imputation_df = (
        df.groupby(["Modified sequence", "Charge"])["MS1 corrected reference intensity"]
        .apply(lambda x: x.mean())
        .reset_index(name="Mean MS1 corrected reference intensity")
    )

    df = pd.merge(
        left=df, right=imputation_df, on=["Modified sequence", "Charge"], how="left"
    )

    has_reference = ~np.isnan(ref_channels.mean(axis=1))
    df.loc[has_reference, "MS1"] = (
        df["Mean MS1 corrected reference intensity"]
        * tmt_channels.sum(axis=1)
        / ref_channels.mean(axis=1)
    )
    df.loc[~has_reference, "MS1"] = np.nan

    return df


def _impute_ms1_intensity(
    df: pd.DataFrame, ref_channel_df: pd.DataFrame
) -> pd.DataFrame:
    """Impute an MS1 value for scans without MS1 but with reference channel values.
    For each modified peptide in each batch if MS1 missing check how many other batches have MS1 and reference channels.
    For those, compute the mean intensity of the reference channels after MS1 redistribution.
    :param df: dataframe of batch intensities (11 reporter channels, intensities are not log transformed and missing values as NaNs
                                               MS1 column, intensities are not log transformed and missing values as zeroes)
    """
    logger.info("Imputing missing MS1 intensities using reference channels")
    df["MS1"] = df["MS1"].replace(0, np.nan)

    tmt_channels = get_tmt_channels(df)
    ref_channels = get_ref_channels(df, ref_channel_df)

    df["MS1 corrected reference intensity"] = (
        ref_channels.mean(axis=1) / tmt_channels.sum(axis=1) * df["MS1"]
    )

    # TODO: switch to ['Modified sequence', 'Charge']
    imputation_df = (
        df.groupby("Modified sequence")["MS1 corrected reference intensity"]
        .apply(lambda x: x.mean())
        .reset_index(name="Mean MS1 corrected reference intensity")
    )

    df = pd.merge(left=df, right=imputation_df, on="Modified sequence", how="left")

    missing_ms1 = np.isnan(df["MS1"])
    has_reference = ~np.isnan(ref_channels.mean(axis=1))
    df.loc[missing_ms1 & has_reference, "MS1"] = (
        df["Mean MS1 corrected reference intensity"]
        * tmt_channels.sum(axis=1)
        / ref_channels.mean(axis=1)
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


def _redistribute_ms1_intensity(df: pd.DataFrame) -> pd.DataFrame:
    """Multiply MS1 with share of intensity in the channel relative to the summed intensity.
    :param df: dataframe of batch intensities (11 reporter channels, intensities are not log transformed and missing values are NaNs
                                               MS1 column, intensities are not log transformed and missing values are zeroes)
    :return:
    """
    logger.info("Distributing MS1 intensity over TMT channels")
    has_ms1 = ~np.isnan(df["MS1"])
    tmt_channels = get_tmt_channels(df)
    df.loc[has_ms1, tmt_channels.columns] = tmt_channels.loc[has_ms1, :].multiply(
        df.loc[has_ms1, "MS1"] / tmt_channels.loc[has_ms1, :].sum(axis=1), axis="index"
    )

    return df
