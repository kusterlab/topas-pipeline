import os
import re
import sys
from pathlib import Path
import logging
import requests
import json
from typing import List
from enum import Enum

import pandas as pd
import numpy as np

from . import config

PATIENT_PREFIX = "pat_"
REF_CHANNEL_PREFIX = "ref_"


# remember to update the corresponding constants in topas-portal/flask-backend/topas_portal/utils.py
class DataType(str, Enum):
    FULL_PROTEOME = "protein"
    FULL_PROTEOME_ANNOTATED = "protein_annotated"
    FULL_PROTEOME_NUM_PEPTIDES = "num_peptides"
    PHOSPHO_PROTEOME = "psite"
    FP_PP = "FP_PP"
    PHOSPHO_PROTEOME_ANNOTATED = "psite_annotated"
    PHOSPHO_SCORE = "phospho_score"
    PHOSPHO_SCORE_PSITE = "phospho_psite"
    KINASE_SCORE = "kinase"
    KINASE_SUBSTRATE = "kinase_substrate"
    TOPAS_KINASE_SCORE = "topas_kinase"
    TOPAS_KINASE_SUBSTRATE = "topas_kinase_substrate"
    TOPAS_PHOSPHO_SCORE = "topas_phospho"
    TOPAS_PHOSPHO_SCORE_PSITE = (
        "topas_phospho_psite"  # p-sites making up a phosphoprotein score
    )
    TOPAS_PROTEIN = "topas_expression"
    TOPAS_RTK_SCORE = "topas_rtk"
    TOPAS_CK_SCORE = "topas_ck"
    TOPAS_SUBSCORE = "topas_subscore"
    BIOMARKER = "biomarker"
    REPORT_SUMMARY = "report_summary"

    PATIENT_METADATA = "patients_df"
    SAMPLE_ANNOTATION = "sample_annotation_df"
    SEARCH_QC = "search_qc"

    TRANSCRIPTOMICS = "fpkm"
    GENOMICS = "genomics"


INDEX_COLS: dict[DataType, list[str]] = {
    DataType.FULL_PROTEOME: ["Gene names"],
    DataType.FULL_PROTEOME_ANNOTATED: ["Gene names"],
    DataType.PHOSPHO_SCORE: ["Gene names"],
    DataType.PHOSPHO_PROTEOME: ["Gene names", "Modified sequence", "Proteins"],
    DataType.PHOSPHO_PROTEOME_ANNOTATED: [
        "Gene names",
        "Modified sequence group",
        "Proteins",
    ],
}


def init_file_logger(results_folder, log_file_name):
    module_name = ".".join(__name__.split(".")[:-1])
    file_logger = logging.FileHandler(results_folder / Path(log_file_name))
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    file_logger.setFormatter(formatter)
    logging.getLogger(module_name).addHandler(file_logger)


def send_slack_message(message: str, results_folder: str, slack_config: config.Slack):
    if slack_config.webhook_url == "":
        return

    message = f"Results folder: {Path(results_folder).name}\n{message}"

    slack_data = {
        "username": "TOPAS pipeline",
        "icon_emoji": ":gem:",
        "channel": slack_config.channel,
        "attachments": [
            {"fields": [{"title": "New Incoming Message", "value": message}]}
        ],
    }
    response = requests.post(
        slack_config.webhook_url,
        data=json.dumps(slack_data),
        headers={
            "Content-Type": "application/json",
            "Content-Length": str(sys.getsizeof(slack_data)),
        },
    )
    if response.status_code != 200:
        raise ValueError(
            "Request to slack returned an error %s, the response is:\n%s"
            % (response.status_code, response.text)
        )


def get_index_cols(data_type: str) -> List[str]:
    index_cols = ["Gene names"]
    if data_type.startswith("pp"):
        index_cols = ["Gene names", "Modified sequence", "Proteins"]
    return index_cols


def short_phospho_notation(modified_sequence_column: pd.Series) -> pd.Series:
    return modified_sequence_column.str.replace(
        r"([STY])\(Phospho \(STY\)\)", lambda pat: f"p{pat.group(1)}", regex=True
    )


def split_str_list(samples_for_report):  # add typehint
    samples_for_report = re.split(r"[;|,]", samples_for_report)
    samples_for_report = [sample.strip() for sample in samples_for_report]
    return samples_for_report


def whitespace_remover(df):
    for col in df.columns:
        if df[col].dtype == "object":
            # applying strip function on column
            df[col] = df[col].astype(str).map(str.strip)
            df[col] = df[col].replace(r"^nan$", np.nan, regex=True)
        else:
            pass
    return df


def get_tmt_channels(df: pd.DataFrame) -> pd.DataFrame:
    return df.filter(regex=r"^Reporter intensity corrected \d{1,2}")


def get_ref_channels(df: pd.DataFrame, ref_channel_df: pd.DataFrame):
    # TODO: add support for different reference channels in each batch
    if (
        ref_channel_df["Batch Name"].nunique() * ref_channel_df["TMT Channel"].nunique()
    ) != len(ref_channel_df[["Batch Name", "TMT Channel"]].drop_duplicates()):
        raise ValueError(
            "Datasets where reference channels differ between batches are not yet supported."
        )

    return df.filter(
        [
            f"Reporter intensity corrected {channel}"
            for channel in ref_channel_df["TMT Channel"].unique()
        ]
    )


def keep_only_sample_columns(df: pd.DataFrame) -> pd.DataFrame:
    return df.filter(regex=rf"^({PATIENT_PREFIX})|(Reporter intensity corrected)")


def keep_only_sample_and_ref_columns(df: pd.DataFrame) -> pd.DataFrame:
    return df.filter(regex=rf"(^{PATIENT_PREFIX})|(^{REF_CHANNEL_PREFIX})")


def add_patient_prefix(patient_list: list[str]):
    return [
        PATIENT_PREFIX + x if not x.startswith(REF_CHANNEL_PREFIX) else x
        for x in patient_list
    ]


def explode_on_separated_string(df: pd.DataFrame, index_col: str, sep: str = ";"):
    index_col_exploded = f"{index_col}_exploded"
    df[index_col_exploded] = df[index_col].str.split(sep)
    return df.explode(index_col_exploded), index_col_exploded


def csv_unique(s: str) -> str:
    return ";".join(sorted(set([x for x in s.split(";") if len(x) > 0])))


def merge_by_delimited_field(
    df: pd.DataFrame,
    other_df: pd.DataFrame,
    field_name: str,
    delimiter: str = ";",
    agg_func=None,
    inplace: bool = False,
) -> pd.DataFrame:
    """
    Merges two dataframe by a field which contains a delimited field in the left dataframe
    """
    if not agg_func:
        agg_func = lambda x: delimiter.join(x.dropna())

    df_with_row_idx = df.assign(row_id=range(len(df)))
    if field_name not in df.columns and field_name in df.index.names:
        df_with_row_idx = df_with_row_idx.reset_index()

    key_df = df_with_row_idx[[field_name, "row_id"]]
    merged_df = (
        key_df.assign(exploded_field=key_df[field_name].str.split(delimiter))
        .explode("exploded_field")
        .drop(columns=field_name)
        .merge(other_df, left_on="exploded_field", right_on=field_name, how="left")
        .drop(columns=[field_name, "exploded_field"])
        .groupby("row_id")
        .agg(agg_func)
        .reset_index()
    )
    if inplace:
        merged_df = df_with_row_idx[["row_id"]].merge(
            merged_df, on="row_id", how="left"
        )
        merged_df = merged_df.drop(columns="row_id")
        merged_df.index = df.index
        df.loc[:, merged_df.columns] = merged_df
    else:
        merged_df = df_with_row_idx.merge(merged_df, on="row_id", how="left")
        merged_df = merged_df.drop(columns="row_id")
        return merged_df


def validate_file_access(func):
    """
    Decorator to check if the file exists and can be opened.

    Usage:

        @validate_file_access
        def your_function(path, ...):
            <code>
    """

    def wrapper(path: str, *args, **kwargs):
        if not os.path.exists(path):
            raise FileNotFoundError(f"{path} does not exist")

        try:
            df = func(path, *args, **kwargs)
            return df
        except PermissionError as err:
            raise PermissionError(
                f"{type(err).__name__}: {err} in reading the file {path}. Check if it is opened somewhere"
            )

    return wrapper
