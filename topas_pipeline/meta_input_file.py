from pathlib import Path
from typing import List, Optional

import pandas as pd


def get_meta_input_file_path(results_folder: str, data_type: str) -> Path:
    return Path(results_folder) / Path(f"meta_input_file_{data_type.upper()}.tsv")


def write_meta_input_file(
    meta_input_file: Path,
    mq_txt_folders: List[str],
    simsi_raw_folders: List[str],
    tmt_correction_files: Optional[List[str]],
):
    meta_input_file_dict = {
        "mq_txt_folder": mq_txt_folders,
        "raw_folder": simsi_raw_folders,
    }
    if tmt_correction_files is not None:
        meta_input_file_dict["tmt_correction_file"] = tmt_correction_files
    meta_input_file_df = pd.DataFrame(meta_input_file_dict)

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


def check_metafiles(results_folder):
    results_folder = Path(results_folder)
    meta_fp = read_meta_input_file(results_folder / "meta_input_file_FP.tsv")
    meta_pp = read_meta_input_file(results_folder / "meta_input_file_PP.tsv")

    # ADD potential for CL in between batch and digits
    fp_batches = meta_fp["mq_txt_folder"].str.extract(r"Batch(\d+)")
    pp_batches = meta_pp["mq_txt_folder"].str.extract(r"Batch(\d+)")

    fp_batches_set = set(fp_batches[0].dropna())
    pp_batches_set = set(pp_batches[0].dropna())

    if fp_batches_set != pp_batches_set:
        raise ValueError(
            "The batches for FP and PP differ (meta_input_file content mismatch)."
        )
