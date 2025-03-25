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


def clinical_process(*args, **kwargs) -> None:
    data_types = kwargs.pop("data_types")
    for data_type in data_types:
        clinical_process_data_type(*args, **kwargs, data_type=data_type)


def clinical_process_data_type(
    results_folder: Union[str, Path],
    debug: bool,
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

    # check if this step has already been done
    if os.path.exists(os.path.join(results_folder, f"annot_{data_type}.csv")):
        logger.info(
            f"Clinical processing {data_type} skipped - found files already preprocessed"
        )
        return
    dfs = dict()
    if data_type == "fp":
        index_col = "Gene names"
        keep_default_na = True
    else:
        index_col = "Modified sequence"
        keep_default_na = False
    dfs[data_type] = pd.read_csv(
        os.path.join(results_folder, f"preprocessed_{data_type}.csv"),
        index_col=index_col,
        keep_default_na=keep_default_na,
    )

    data_type_with_ref = f"{data_type}_with_ref"
    dfs[data_type_with_ref] = pd.read_csv(
        os.path.join(results_folder, f"preprocessed_{data_type_with_ref}.csv"),
        index_col=index_col,
        keep_default_na=keep_default_na,
    )

    if data_type == "pp":
        logger.info("Annotating phospho sites")
        dfs[data_type] = clinical_tools.phospho_annot(
            dfs[data_type],
            clinic_proc_config,
        )

        logger.info("Adding PhosphoSitePlus URLs")
        # TODO: add PSP URLs again
        dfs[data_type] = clinical_tools.add_psp_urls(dfs[data_type])

        if debug:
            dfs[data_type_with_ref] = clinical_tools.phospho_annot(
                dfs[data_type_with_ref],
                clinic_proc_config,
            )
            dfs[data_type_with_ref] = clinical_tools.add_psp_urls(
                dfs[data_type_with_ref]
            )

    annot_levels = list(TOPAS_SCORE_COLUMNS.keys()) + list(
        TOPAS_SUBSCORE_COLUMNS.keys()
    )

    for data_type in dfs.keys():
        for annot_type in annot_levels:
            dfs[data_type] = clinical_tools.prot_clinical_annotation(
                dfs[data_type],
                clinic_proc_config.prot_baskets,
                data_type=data_type,
                annot_type=annot_type,
            )
        dfs[data_type].to_csv(os.path.join(results_folder, f"annot_{data_type}.csv"), float_format="%.6g")


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


def read_annotation_files(
    results_folder: Union[str, Path],
    debug: bool,
    data_type: str,
    post_process_topas_columns: bool = True,
):
    index_cols = utils.get_index_cols(data_type)

    annot = pd.read_csv(
        os.path.join(results_folder, f"annot_{data_type}.csv"), index_col=index_cols
    )

    annot_ref = None
    if debug:
        annot_ref = pd.read_csv(
            os.path.join(results_folder, f"annot_{data_type}_with_ref.csv"),
            index_col=index_cols,
        )
        if post_process_topas_columns:
            for key in TOPAS_SUBSCORE_COLUMNS.keys():
                annot_ref[key] = annot_ref.apply(merge_topas_score_and_subscore_names, axis=1)
            for key in TOPAS_SCORE_COLUMNS.keys():
                annot_ref[key] = annot_ref[key].apply(get_unique_topas_names)

    if post_process_topas_columns:
        # Get unique TOPAS names and add main TOPAS name to TOPAS subscore level
        for key in TOPAS_SUBSCORE_COLUMNS.keys():
            annot[key] = annot.apply(merge_topas_score_and_subscore_names, axis=1)
        for key in TOPAS_SCORE_COLUMNS.keys():
            annot[key] = annot[key].apply(get_unique_topas_names)
    return annot, annot_ref


"""
python3 -m topas_pipeline.clinical_process -c config_patients.json
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

    clinical_process(
        results_folder=configs.results_folder,
        debug=configs.preprocessing.debug,
        clinic_proc_config=configs.clinic_proc,
        data_types=configs.data_types,
    )
