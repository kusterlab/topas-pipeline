from typing import Dict, List

import pandas as pd

from .. import identification_metadata
from .. import utils


def rename_columns_with_sample_ids(
    df: pd.DataFrame,
    channel_to_sample_id_dict: Dict[str, str],
    index_cols: List[str],
) -> pd.DataFrame:
    """
    Transform column names of the format 'Reporter intensity corrected <TMT_channel> <batch_name>' to the sample names
    """
    df = df.set_index(index_cols)
    tmt_channels = list(channel_to_sample_id_dict.keys())

    # keep channels that are in the dataframe but missing in the sample annotation file, marked with a prefix
    channel_to_sample_id_dict |= get_missing_channels_dict(df, tmt_channels)
    tmt_channels = list(channel_to_sample_id_dict.keys())

    metadata_cols = list(map(identification_metadata.as_metadata_columns, tmt_channels))

    keep_cols = tmt_channels + metadata_cols
    df = df.filter(items=keep_cols, axis=1)

    # build dictionary to also rename the metadata columns with the sample ids
    metadata_cols_with_sample_id = map(
        lambda x: f"Identification metadata {x}", channel_to_sample_id_dict.values()
    )
    metadata_to_sample_id_dict = dict(zip(metadata_cols, metadata_cols_with_sample_id))

    # add prefix to patient intensity columns
    channel_to_sample_id_dict_final = {
        k: utils.PATIENT_PREFIX + v
        for k, v in channel_to_sample_id_dict.items()
        if not v.startswith(utils.REF_CHANNEL_PREFIX)
        and not v.startswith(utils.OTHER_CHANNEL_PREFIX)
    }
    # keep prefix for intensity columns for reference samples and samples missing from the sample annotation file
    channel_to_sample_id_dict_final |= {
        k: v
        for k, v in channel_to_sample_id_dict.items()
        if v.startswith(utils.REF_CHANNEL_PREFIX)
        or v.startswith(utils.OTHER_CHANNEL_PREFIX)
    }

    rename_dict = {**channel_to_sample_id_dict_final, **metadata_to_sample_id_dict}
    df = df.rename(columns=rename_dict)
    df = df.reset_index()
    return df


def get_missing_channels_dict(df: pd.DataFrame, tmt_channels: list[str]):
    intensity_columns: pd.Index = df.filter(
        regex="^Reporter intensity corrected"
    ).columns
    missing_channels_df: pd.DataFrame = intensity_columns.difference(
        tmt_channels
    ).to_frame(name="Sample name")
    missing_channels_df["batch"] = (
        missing_channels_df["Sample name"].str.split("_Batch").str[-1]
    )
    missing_channels_df["channel"] = (
        missing_channels_df["Sample name"].str.split(" ").str[-2]
    )
    missing_channels_df["Sample name"] = (
        utils.OTHER_CHANNEL_PREFIX
        + "channel_"
        + missing_channels_df["channel"]
        + "_batch"
        + missing_channels_df["batch"]
    )
    return missing_channels_df["Sample name"].to_dict()
