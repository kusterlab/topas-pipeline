import os
import sys
import logging

import pandas as pd
from pathlib import Path
from typing import Union

from . import config
from . import utils
from .topas import protein_phosphorylation
from .annotation import proteins_of_interest as poi
from .annotation import phosphosite
from .preprocess import phospho_grouping

# hacky way to get the package logger instead of just __main__ when running as python -m topas_pipeline.simsi ...
logger = logging.getLogger(__package__ + "." + __file__)


def add_clinical_annotations(*args, **kwargs) -> None:
    data_types: list[str] = kwargs.pop("data_types")
    if "pp" in data_types:
        data_types.append("phospho_score")

    for data_type in data_types:
        add_clinical_annotations_data_type(*args, **kwargs, data_type=data_type)


def add_clinical_annotations_data_type(
    results_folder: Union[str, Path],
    clinic_proc_config: config.ClinicProc,
    data_type: str,
    overwrite: bool = False,
) -> None:
    """
    Opens preprocessed data, annotates phospho and TOPAS scores

    :param results_folder: path to results folder to which preprocessed data can be found and clinically processed data saved
    :param debug:
    :param clinic_proc_config: paths to file used for annotating to phospho and TOPAS scores
    :param data_type: 'fp' for full proteome, 'pp' for phospho proteome
    """
    if os.path.exists(os.path.join(results_folder, f"annot_{data_type}.csv")):
        if not overwrite:
            logger.info(
                f"Clinical processing {data_type} skipped - found files already preprocessed"
            )
            return
        logger.info(f"Found existing results but overwrite flag was set.")

    if data_type == "pp":
        preprocessed_df = phospho_grouping.read_cohort_intensities_df(
            os.path.join(results_folder, "preprocessed_pp.csv"),
            keep_identification_metadata_columns=True,
        )
        annot_df = build_index_annotation_df(preprocessed_df)

        annot_df = phosphosite.add_ck_substrate_annotations(
            annot_df,
            results_folder,
            clinic_proc_config.topas_kinase_substrate_file,
        )

        # ~30 minutes for 3000 samples
        annot_df = phosphosite.add_phospho_annotations(
            annot_df,
            clinic_proc_config,
        )

        annot_df = phosphosite.combine_topas_annotations(annot_df)

        annot_df = merge_intensities_and_annotation_dfs(annot_df, preprocessed_df)
    elif data_type == "phospho_score":
        annot_df = protein_phosphorylation.load_protein_phosphorylation(
            Path(results_folder), protein_results_folder="topas_scores"
        )
    elif data_type == "fp":
        annot_df = pd.read_csv(
            os.path.join(results_folder, f"preprocessed_fp.csv"),
            index_col="Gene names",
        )
    else:
        raise ValueError(f"Unknown data type {data_type}")

    if len(clinic_proc_config.proteins_of_interest_file) > 0:
        logger.info(f"Annotating proteins of interest for {data_type}")
        poi_annotation_df = poi.load_poi_annotation_df(
            clinic_proc_config.proteins_of_interest_file
        )
        poi.merge_with_poi_annotations_inplace(annot_df, poi_annotation_df)

    logger.info(f"Writing annotated table for {data_type}")
    annot_df.to_csv(
        os.path.join(results_folder, f"annot_{data_type}.csv"), float_format="%.6g"
    )


def build_index_annotation_df(preprocessed_df: pd.DataFrame) -> pd.DataFrame:
    keep_cols = phospho_grouping.INDEX_COLS.copy()
    keep_cols.remove("Modified sequence group")

    annot_df = pd.DataFrame(index=preprocessed_df.index)
    annot_df = annot_df.reset_index(keep_cols)

    # annot_df["Modified sequence"] = annot_df.index.str.split(";")
    # annot_df = annot_df.explode("Modified sequence")
    annot_df["Modified sequence"] = annot_df["Modified sequence representative"]
    return annot_df


def merge_intensities_and_annotation_dfs(
    annot_df: pd.DataFrame, preprocessed_df: pd.DataFrame
) -> pd.DataFrame:
    keep_cols = phospho_grouping.INDEX_COLS.copy()
    keep_cols.remove("Modified sequence group")

    annot_df = (
        annot_df.reset_index()
    )  # make "Modified sequence a regular column to keep it during merging

    annot_df = preprocessed_df.reset_index(keep_cols, drop=True).merge(
        annot_df, on="Modified sequence group"
    )
    annot_df = annot_df.set_index("Modified sequence group")
    return annot_df


def read_annotated_expression_file(
    results_folder: Union[str, Path],
    data_type: str,
) -> pd.DataFrame:
    if data_type == "pp":
        data_type_enum = utils.DataType.PHOSPHO_PROTEOME_ANNOTATED
    elif data_type == "fp":
        data_type_enum = utils.DataType.FULL_PROTEOME_ANNOTATED
    elif data_type == "phospho_score":
        data_type_enum = utils.DataType.PHOSPHO_SCORE
    else:
        raise ValueError(
            f"Unknown data type '{data_type}' for read_annotated_expression_file"
        )
    return pd.read_csv(
        os.path.join(results_folder, f"annot_{data_type}.csv"),
        index_col=utils.INDEX_COLS[data_type_enum],
    )


"""
python3 -m topas_pipeline.clinical_annotation -c config_patients.json
"""
if __name__ == "__main__":
    import argparse

    from . import config

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-c", "--config", required=True, help="Absolute path to configuration file."
    )
    parser.add_argument(
        "-o",
        "--overwrite",
        action="store_true",
        help="Ignore existing results and recompute outputs.",
    )
    args = parser.parse_args(sys.argv[1:])

    configs = config.load(args.config)

    add_clinical_annotations(
        results_folder=configs.results_folder,
        clinic_proc_config=configs.clinic_proc,
        data_types=configs.data_types,
        overwrite=args.overwrite,
    )
