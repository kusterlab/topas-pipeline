from pathlib import Path
from typing import List

import pandas as pd


def get_meta_input_file_path(results_folder: str, data_type: str) -> Path:
    return Path(results_folder) / Path(f"meta_input_file_{data_type.upper()}.tsv")


def write_meta_input_file(
    meta_input_file: Path,
    mq_txt_folders: List[str],
    simsi_raw_folders: List[str],
    tmt_correction_files: List[str],
):
    meta_input_file_df = pd.DataFrame(
        {
            "mq_txt_folder": mq_txt_folders,
            "raw_folder": simsi_raw_folders,
            "tmt_correction_file": tmt_correction_files,
        }
    )

    meta_input_file_df.to_csv(meta_input_file, sep="\t", index=False)


def meta_input_files_equal(meta_input_file: Path, meta_input_file_other: Path) -> bool:
    meta_input_df = read_meta_input_file(meta_input_file)
    meta_input_other_df = read_meta_input_file(meta_input_file_other)
    return meta_input_df.sort_values(by="mq_txt_folder").equals(
        meta_input_other_df.sort_values(by="mq_txt_folder")
    )


def read_meta_input_file(meta_input_file: Path) -> pd.DataFrame:
    return pd.read_csv(
        meta_input_file, sep="\t", usecols=["mq_txt_folder", "raw_folder"]
    )
