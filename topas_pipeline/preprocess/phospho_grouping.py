import os
import sys
import logging
import collections
from pathlib import Path

import pandas as pd
import numpy as np
import psite_annotation as pa

from tqdm import tqdm

tqdm.pandas()

from .. import sample_annotation
from . import sample_mapping
from .. import identification_metadata

# hacky way to get the package logger instead of just __main__ when running as a module
logger = logging.getLogger(__package__ + "." + Path(__file__).stem)


INDEX_COLS = [
    "Modified sequence group",
    "Modified sequence representative",
    "Gene names",
    "Proteins",
]


def aggregate_modified_sequences(results_folder: str) -> pd.DataFrame:
    """Aggregate modified sequence rows with localization within +/-2 amino acids.

    Args:
        results_folder (str): _description_
    """
    results_folder = Path(results_folder)
    if os.path.exists(results_folder / "preprocessed_pp2_agg.csv"):
        logger.info(f"Phospho grouping skipped - found file already processed")
        return
    pp_df = read_preprocessed_pp2(results_folder)

    # clean Mod_seq string
    replace_dict = {
        r"\(Acetyl \(Protein N-term\)\)": "(ac)",
        r"M\(Oxidation \(M\)\)": "M",
        r"p([STY])": r"\1(ph)",
    }

    pp_df["Modified sequence"] = pp_df["Modified sequence"].str[1:-1]
    for old, new in replace_dict.items():
        pp_df["Modified sequence"] = pp_df["Modified sequence"].replace(
            old, new, regex=True
        )

    # replace empty metadata cells (=measured in sample) with ";" to recognize
    # when imputed data is combined with measured values in the aggregation
    # 1.5 minutes for full matrix
    logger.info("Aggregating modified sequence groups")
    pp_df.loc[:, pp_df.filter(like="Identification metadata").columns] = pp_df.loc[
        :, pp_df.filter(like="Identification metadata").columns
    ].fillna(";")

    # 4.5 minutes for full matrix
    agg_pp_df = pa.aggregateModifiedSequenceGroups(
        pp_df,
        experiment_cols=pp_df.filter(like="Reporter intensity corrected").columns,
        agg_cols={"Gene names": "first", "Proteins": "first"}
        | {c: "sum" for c in pp_df.filter(like="Identification metadata").columns},
    )

    # aggregate metadata imputation status. the most common case is a non-aggregated
    # measured value, we convert these into nans to speed up the .map() function
    # 3.5 minutes for full matrix
    agg_pp_df.loc[:, agg_pp_df.filter(like="Identification metadata").columns] = (
        agg_pp_df.loc[:, agg_pp_df.filter(like="Identification metadata").columns]
        .replace(r"^;+$|^$", np.nan, regex=True)
        .map(aggregate_imputations, na_action="ignore")
    )

    # 7.5 minutes for full matrix
    agg_pp_file = results_folder / "preprocessed_pp2_agg.csv"
    logger.info(f"Writing aggregated modified sequence groups to {agg_pp_file}")
    agg_pp_df.to_csv(agg_pp_file, index=False, float_format="%.6g")

    return agg_pp_df


def read_preprocessed_pp2(results_folder: str) -> pd.DataFrame:
    logger.info("Reading in preprocessed_pp2.csv")
    headers = pd.read_csv(results_folder / "preprocessed_pp2.csv", nrows=1)
    dtype_dict = collections.defaultdict(
        lambda: "str"
    )  # 'str' dtype does still create NaNs for empty cells
    dtype_dict |= {
        c: "float32"
        for c in headers.filter(like="Reporter intensity corrected").columns
    }

    pp_df = pd.read_csv(results_folder / "preprocessed_pp2.csv", dtype=dtype_dict)
    return pp_df


def aggregate_imputations(x):
    annotations = set(x[:-1].split(";"))

    if len(annotations) >= 2:
        annotations = summarize_annotations(annotations)
        return annotations
    # else, just one annotation
    return ";".join(map(str, list(annotations))) + ";"


def summarize_annotations(annotations):
    annotations = set(annotations)  # ensure set

    if not has_quan_OOR(annotations):
        # just a fallback logic. Shoult not happen
        if not has_partially_imputed(annotations) and not has_empty(annotations):
            return "imputed;"
        return "partially imputed;"

    # has quan_OOR
    annotations.discard("quan_OOR")

    if not has_imputed(annotations) and not has_partially_imputed(annotations):
        return "partially quan_OOR;"
    elif has_partially_imputed(annotations) or (
        has_imputed(annotations) and has_empty(annotations)
    ):
        return "partially imputed;quan_OOR;"
    elif annotations == {"imputed"}:
        return "imputed;quan_OOR;"
    else:
        raise ValueError(f"Unexpected annotation combination: {annotations}")


def has_quan_OOR(annotations):
    return "quan_OOR" in annotations


def has_imputed(annotations):
    return "imputed" in annotations


def has_partially_imputed(annotations):
    return "partially imputed" in annotations


def has_empty(annotations):
    return "" in annotations


def read_cohort_intensities_df(
    grouped_phospho_file: str,
    sample_annotation_file: str = None,
    skiprows: pd.Series = None,
    keep_identification_metadata_columns: bool = False,
):
    """
    Read in preprocessed_pp.csv or preprocessed_pp2_agg.csv
    """
    logger.info(f"Reading {Path(grouped_phospho_file).name}")
    header = pd.read_csv(grouped_phospho_file, index_col=0, nrows=1)
    intensity_cols = header.filter(like="Reporter intensity corrected").columns.tolist()
    if len(intensity_cols) == 0:
        intensity_cols = header.filter(regex=r"(^pat_)|(^ref_)").columns.tolist()

    identification_metadata_cols = []
    if keep_identification_metadata_columns:
        identification_metadata_cols = header.filter(
            like=identification_metadata.METADATA_COLUMN_PREFIX
        ).columns.tolist()

    dtype_dict = collections.defaultdict(lambda: "string")
    dtype_dict |= {c: "float64" for c in intensity_cols}
    dtype_dict |= {c: "string" for c in identification_metadata_cols}

    intensities_df = pd.read_csv(
        grouped_phospho_file,
        skiprows=skiprows,
        usecols=INDEX_COLS + intensity_cols + identification_metadata_cols,
        dtype=dtype_dict,
    )
    # fill na values in the Gene names and Proteins columns
    intensities_df[INDEX_COLS] = intensities_df[INDEX_COLS].fillna("")

    intensities_df = intensities_df.set_index(INDEX_COLS)

    if sample_annotation_file:
        sample_annotation_df = pd.read_csv(sample_annotation_file)
        channel_to_sample_id_dict = sample_annotation.get_channel_to_sample_id_dict(
            sample_annotation_df,
            remove_qc_failed=True,
            remove_replicates=False,
        )

        intensities_df = sample_mapping.rename_columns_with_sample_ids(
            intensities_df.reset_index(),
            channel_to_sample_id_dict,
            index_cols=INDEX_COLS,
        )
        intensities_df = intensities_df.set_index(INDEX_COLS)

    return intensities_df


"""
python3 -m topas_pipeline.phospho_grouping -c config_patients.json
"""
if __name__ == "__main__":
    import argparse

    from . import config

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-c", "--config", required=True, help="Absolute path to configuration file."
    )
    args = parser.parse_args(sys.argv[1:])

    configs = config.load(args.config)

    aggregate_modified_sequences(
        results_folder=configs.results_folder,
    )
