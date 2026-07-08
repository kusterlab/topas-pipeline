import numpy as np
import pandas as pd
import re
from typing import Union


def keep_protein_name_from_fasta_header(value: Union[str, float]):
    if not pd.isnull(value):
        if "sp" in value:
            return value.split("|")[1]
        else:
            return value
    else:
        return np.nan


def convert_mods_to_mq_format(
    mod_seq_series: pd.Series,
    replacements: dict[str, str],
    add_flanking_underscores: bool = True,
):
    for old, new in replacements.items():
        mod_seq_series = mod_seq_series.str.replace(old, new, regex=False)

    if add_flanking_underscores:
        mod_seq_series = "_" + mod_seq_series + "_"
    return mod_seq_series


def convert_mod_column_to_mq_format(fragpipe_mod_str: str) -> str:
    if "79.9663" in fragpipe_mod_str:
        return "Phospho (STY)"
    return ""


def modify_protein_and_peptide_info(df: pd.DataFrame) -> pd.DataFrame:
    df["Proteins"] = df["Proteins"].apply(keep_protein_name_from_fasta_header)
    df["Leading proteins"] = df["Leading proteins"].apply(
        keep_protein_name_from_fasta_header
    )

    # pretend that the LFQ data is TMT data with a single TMT channel
    df["Reporter intensity corrected 1"] = df["Intensity"]

    # convert phospho modification notation, e.g. S(Phospho (STY)) => pS
    df["Modified sequence"] = df["Modified sequence"].str.replace(
        re.compile(r"([STY])\(Phospho \(STY\)\)"),
        lambda pat: f"p{pat.group(1)}",
        regex=True,
    )

    return df
