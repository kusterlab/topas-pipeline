from typing import Dict, List

import pandas as pd

from .. import identification_metadata


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
    metadata_cols = list(map(identification_metadata.as_metadata_columns, tmt_channels))

    keep_cols = tmt_channels + metadata_cols
    df = df.filter(items=keep_cols, axis=1)

    # build dictionary to also rename the metadata columns with the sample ids
    metadata_cols_with_sample_id = map(
        lambda x: f"Identification metadata {x}", channel_to_sample_id_dict.values()
    )
    metadata_to_sample_id_dict = dict(zip(metadata_cols, metadata_cols_with_sample_id))

    # add "pat_" prefix to patient intensity columns
    channel_to_sample_id_dict_final = {
        k: "pat_" + v
        for k, v in channel_to_sample_id_dict.items()
        if not v.startswith("ref_")
    }
    # keep "ref_" prefix for reference intensity columns
    channel_to_sample_id_dict_final |= {
        k: v for k, v in channel_to_sample_id_dict.items() if v.startswith("ref_")
    }

    rename_dict = {**channel_to_sample_id_dict_final, **metadata_to_sample_id_dict}
    df = df.rename(columns=rename_dict)
    df = df.reset_index()
    return df
