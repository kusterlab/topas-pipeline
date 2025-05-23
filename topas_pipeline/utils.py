import os
import re
import sys
from pathlib import Path
import logging
import requests
import time
import json
from typing import List

import pandas as pd
import numpy as np

from . import config


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
    return df.filter(regex=r"^(pat_)|(Reporter intensity corrected)")


def explode_on_separated_string(df: pd.DataFrame, index_col: str, sep: str = ";"):
    index_col_exploded = f"{index_col}_exploded"
    df[index_col_exploded] = df[index_col].str.split(sep)
    return df.explode(index_col_exploded), index_col_exploded


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
