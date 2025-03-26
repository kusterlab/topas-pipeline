import os
import sys
import logging

import numpy as np
import pandas as pd
from pathlib import Path
from typing import List, Union

from . import config
from . import clinical_tools
from . import utils

logger = logging.getLogger(__name__)

TOPAS_SCORE_COLUMNS = {
    "TOPAS_score": "TOPAS annot",
    "POI_category": "POI category",
}
TOPAS_SUBSCORE_COLUMNS = {
    "TOPAS_subscore": "TOPAS sublevel annot",
}

# # Old annotation file format column names
# TOPAS_SCORE_COLUMNS = {
#     "basket": "TOPAS annot",
# }
# TOPAS_SUBSCORE_COLUMNS = {
#     "sub_basket": "TOPAS sublevel annot",
# }


def add_clinical_annotations(*args, **kwargs) -> None:
    data_types = kwargs.pop("data_types")
    debug = kwargs.pop("debug")
    if debug:
        data_types += [data_type + "_with_ref" for data_type in data_types]

    for data_type in data_types:
        add_clinical_annotations_data_type(*args, **kwargs, data_type=data_type)


def add_clinical_annotations_data_type(
    results_folder: Union[str, Path],
    clinic_proc_config: config.ClinicProc,
    data_type: str,
) -> None:
    """
    Opens preprocessed data, annotates phospho and TOPAS scores

    :param results_folder: path to results folder to which preprocessed data can be found and clinically processed data saved
    :param debug:
    :param clinic_proc_config: paths to file used for annotating to phospho and TOPAS scpres
    :param data_type: 'fp' for full proteome, 'pp' for phospho proteome
    """
    # TODO: check if data files with ref if so use these as with_ref otherwise use normal as with_ref

    if os.path.exists(os.path.join(results_folder, f"annot_{data_type}.csv")):
        logger.info(
            f"Clinical processing {data_type} skipped - found files already preprocessed"
        )
        return

    if data_type.startswith("fp"):
        index_col = "Gene names"
        keep_default_na = True
    else:
        index_col = "Modified sequence"
        keep_default_na = False
    preprocessed_df = pd.read_csv(
        os.path.join(results_folder, f"preprocessed_{data_type}.csv"),
        index_col=index_col,
        keep_default_na=keep_default_na,
    )

    if data_type == "pp":
        logger.info("Annotating phospho sites")
        preprocessed_df = clinical_tools.phospho_annot(
            preprocessed_df,
            clinic_proc_config,
        )

        logger.info("Adding PhosphoSitePlus URLs")
        preprocessed_df = clinical_tools.add_psp_urls(preprocessed_df)

    annot_levels = list(TOPAS_SCORE_COLUMNS.keys()) + list(
        TOPAS_SUBSCORE_COLUMNS.keys()
    )

    for annot_type in annot_levels:
        preprocessed_df = clinical_tools.prot_clinical_annotation(
            preprocessed_df,
            clinic_proc_config.prot_baskets,
            data_type=data_type,
            annot_type=annot_type,
        )
    preprocessed_df.to_csv(
        os.path.join(results_folder, f"annot_{data_type}.csv"), float_format="%.6g"
    )


def merge_topas_score_and_subscore_names(row: pd.Series) -> str:
    topas_subscore_names = row[list(TOPAS_SUBSCORE_COLUMNS.keys())[0]]
    if not pd.isnull(row[list(TOPAS_SUBSCORE_COLUMNS.keys())[0]]):
        topas_subscore_list = row[list(TOPAS_SUBSCORE_COLUMNS.keys())[0]].split(";")
        topas_score_list = row[list(TOPAS_SCORE_COLUMNS.keys())[0]].split(";")
        topas_subscore_names = [
            (
                topas_score_list[i] + " - " + topas_subscore_list[i]
                if len(topas_score_list[i]) > 0
                else ""
            )
            for i in range(len(topas_subscore_list))
        ]
        topas_subscore_names = get_unique_topas_names(topas_subscore_names)
    return topas_subscore_names


def get_unique_topas_names(topas_names: Union[List, str, float]) -> str:
    if type(topas_names) != list and type(topas_names) != float:
        topas_names = topas_names.split(";")
    if type(topas_names) != float:
        topas_names = ";".join(np.unique(np.array(topas_names)))
    return topas_names


def read_annotated_expression_file(
    results_folder: Union[str, Path],
    data_type: str,
) -> pd.DataFrame:
    return pd.read_csv(
        os.path.join(results_folder, f"annot_{data_type}.csv"),
        index_col=utils.get_index_cols(data_type),
    )


def post_process_topas_columns(annot: pd.DataFrame) -> pd.DataFrame:
    # Get unique TOPAS names and add main TOPAS name to TOPAS subscore level
    for key in TOPAS_SUBSCORE_COLUMNS.keys():
        annot[key] = annot.apply(merge_topas_score_and_subscore_names, axis=1)
    for key in TOPAS_SCORE_COLUMNS.keys():
        annot[key] = annot[key].apply(get_unique_topas_names)
    return annot


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
    args = parser.parse_args(sys.argv[1:])

    configs = config.load(args.config)

    add_clinical_annotations(
        results_folder=configs.results_folder,
        debug=configs.preprocessing.debug,
        clinic_proc_config=configs.clinic_proc,
        data_types=configs.data_types,
    )
