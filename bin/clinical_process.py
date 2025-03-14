import os
import sys
import json
import logging

import numpy as np
import pandas as pd
from pathlib import Path
from typing import List, Union

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
    data_types = kwargs.pop('data_types')
    for data_type in data_types:
        clinical_process_data_type(*args, **kwargs, data_type=data_type)


def clinical_process_data_type(results_folder: Union[str, Path],
                               debug: bool,
                               prot_baskets: Union[str, Path],
                               extra_kinase_annot: Union[str, Path],
                               pspFastaFile: Union[str, Path],
                               pspKinaseSubstrateFile: Union[str, Path],
                               pspAnnotationFile: Union[str, Path],
                               pspRegulatoryFile: Union[str, Path],
                               data_type: str) -> None:
    """
    Opens preprocessed data, annotates phospho and clinical baskets

    :param results_folder: location of results folder to which preprocessed data can be found and clinically processed data saved
    :param debug:
    :param prot_baskets: path to file used for annotating to clinical baskets/pathways
    :param extra_kinase_annot: path to file with extra kinase annotations (optional)
    :param pspFastaFile: file used for adding peptide and psite positions
    :param pspKinaseSubstrateFile: file used for adding kinase substrate annotation
    :param pspAnnotationFile: file used for adding annotations from Phosphosite Plus
    :param pspRegulatoryFile: file used for adding regulatory information
    :param data_type:
    """
    # TODO: check if data files with ref if so use these as with_ref otherwise use normal as with_ref

    # check if this step has already been done
    if os.path.exists(os.path.join(results_folder, f'annot_{data_type}.csv')):
        logger.info(f'Clinical processing {data_type} skipped - found files already preprocessed')
        return
    dfs = dict()
    if data_type == 'fp':
        index_col = 'Gene names'
        keep_default_na = True
    else:
        index_col = 'Modified sequence'
        keep_default_na = False
    dfs[data_type] = pd.read_csv(os.path.join(results_folder, f'preprocessed_{data_type}.csv'),
                                 index_col=index_col, keep_default_na=keep_default_na)

    data_type_with_ref = f'{data_type}_with_ref'
    dfs[data_type_with_ref] = pd.read_csv(os.path.join(results_folder, f'preprocessed_{data_type_with_ref}.csv'),
                                          index_col=index_col, keep_default_na=keep_default_na)

    if data_type == 'pp':
        logger.info('Annotating phospho sites')
        dfs[data_type] = clinical_tools.phospho_annot(dfs[data_type], extra_kinase_annot, pspFastaFile, pspKinaseSubstrateFile,
                                     pspAnnotationFile, pspRegulatoryFile)

        logger.info('Adding PhosphoSitePlus URLs')
        # TODO: add PSP URLs again
        dfs[data_type] = clinical_tools.add_psp_urls(dfs[data_type])

        if debug:
            dfs[data_type_with_ref] = clinical_tools.phospho_annot(dfs[data_type_with_ref], extra_kinase_annot, pspFastaFile, pspKinaseSubstrateFile,
                                         pspAnnotationFile, pspRegulatoryFile)
            dfs[data_type_with_ref] = clinical_tools.add_psp_urls(dfs[data_type_with_ref])

    annot_levels = list(TOPAS_SCORE_COLUMNS.keys()) + list(TOPAS_SUBSCORE_COLUMNS.keys())
    
    for data_type in dfs.keys():
        for annot_type in annot_levels:
            dfs[data_type] = clinical_tools.prot_clinical_annotation(dfs[data_type], prot_baskets,
                                                                            data_type=data_type,
                                                                            annot_type=annot_type)
        dfs[data_type].to_csv(os.path.join(results_folder, f'annot_{data_type}.csv'))


def merge_baskets_with_subbaskets(row: pd.Series) -> str:
    subbasket = row[list(TOPAS_SUBSCORE_COLUMNS.keys())[0]]
    if not pd.isnull(row[list(TOPAS_SUBSCORE_COLUMNS.keys())[0]]):
        subbasket_list = row[list(TOPAS_SUBSCORE_COLUMNS.keys())[0]].split(';')
        basket_list = row[list(TOPAS_SCORE_COLUMNS.keys())[0]].split(';')
        subbasket = [basket_list[i] + " - " + subbasket_list[i] if len(basket_list[i]) > 0 else '' for i in
                     range(len(subbasket_list))]
        subbasket = get_unique_baskets(subbasket)
    return subbasket


def get_unique_baskets(baskets: Union[List, str, float]) -> str:
    if type(baskets) != list and type(baskets) != float:
        baskets = baskets.split(';')
    if type(baskets) != float:
        baskets = ";".join(np.unique(np.array(baskets)))
    return baskets


def read_annotation_files(results_folder: Union[str, Path],
                          debug: bool,
                          data_type: str,
                          post_process_basket_columns: bool=True):
    index_cols = utils.get_index_cols(data_type)

    annot = pd.read_csv(os.path.join(results_folder, f'annot_{data_type}.csv'), index_col=index_cols)

    annot_ref = None
    if debug:
        annot_ref = pd.read_csv(os.path.join(results_folder, f'annot_{data_type}_with_ref.csv'), index_col=index_cols)
        if post_process_basket_columns:
            for key in TOPAS_SUBSCORE_COLUMNS.keys():
                annot_ref[key] = annot_ref.apply(merge_baskets_with_subbaskets, axis=1)
            for key in TOPAS_SCORE_COLUMNS.keys():
                annot_ref[key] = annot_ref[key].apply(get_unique_baskets)
    
    if post_process_basket_columns:
        # Get unique baskets and add main basket name to subbasket level
        for key in TOPAS_SUBSCORE_COLUMNS.keys():
            annot[key] = annot.apply(merge_baskets_with_subbaskets, axis=1)
        for key in TOPAS_SCORE_COLUMNS.keys():
            annot[key] = annot[key].apply(get_unique_baskets)
    return annot, annot_ref


"""
python3 -m bin.clinical_process -c config_patients.json
"""
if __name__ == '__main__':
    import argparse

    from . import config

    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", required=True,
                        help="Absolute path to configuration file.")
    args = parser.parse_args(sys.argv[1:])

    configs = config.load(args.config)

    clinical_process(configs["results_folder"], configs["preprocessing"]["debug"], **configs["clinic_proc"],
                     data_types=configs["data_types"])
