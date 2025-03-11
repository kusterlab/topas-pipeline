import sys
import json
import logging
from itertools import compress
from pathlib import Path
from typing import List, Dict, Tuple, Union

import pandas as pd
import numpy as np

import psite_annotation as pa
from . import utils

logger = logging.getLogger(__name__)


def phospho_annot(df: pd.DataFrame,
                  extra_kinase_annot: Union[str, Path, None, int],
                  pspFastaFile: Union[str, Path],
                  pspKinaseSubstrateFile: Union[str, Path],
                  pspAnnotationFile: Union[str, Path],
                  pspRegulatoryFile: Union[str, Path]) -> pd.DataFrame:
    """
    Phospho-site annotation of experimental data using in-house developed tool (MT) based mainly on Phosphosite Plus but also in vitro
    experiments from the lab of Ishihama.

    :param df: dataframe with measured peptide intensities
    :param pspFastaFile: file used for adding peptide and psite positions
    :param pspKinaseSubstrateFile: file used for adding kinase substrate annotation
    :param pspAnnotationFile: file used for adding annotations from Phosphosite Plus
    :param pspRegulatoryFile: file used for adding regulatory information
    """
    logger.info('Phosphosite annotation')

    logger.info(f'Phospho data before: {df.shape}')
    df = df.reset_index()
    df = pa.addPeptideAndPsitePositions(df, pspFastaFile, pspInput=True, returnAllPotentialSites=False)
    df = pa.addPSPKinaseSubstrateAnnotations(df, pspKinaseSubstrateFile, gene_name=True)
    df = pa.addPSPAnnotations(df, pspAnnotationFile)
    df = pa.addPSPRegulatoryAnnotations(df, pspRegulatoryFile)
    df['PSP_LT_LIT'] = df['PSP_LT_LIT'].apply(lambda x: max(x.split(';')))
    df['PSP_MS_LIT'] = df['PSP_MS_LIT'].apply(lambda x: max(x.split(';')))
    df['PSP_MS_CST'] = df['PSP_MS_CST'].apply(lambda x: max(x.split(';')))
    df.rename(columns={'Site positions': 'Site positions identified (MQ)'}, inplace=True)
    df = pa.addPeptideAndPsitePositions(df, pspFastaFile, pspInput=True, returnAllPotentialSites=True)
    logger.info(f'Phospho data after: {df.shape}')
    df = df.set_index('Modified sequence', drop=True)

    # Add extra kinase annotations
    if isinstance(extra_kinase_annot, str):
        df = add_extra_kinase_annotations(df, extra_kinase_annot)
    
    return df


def add_extra_kinase_annotations(df: pd.DataFrame, extra_kinase_annot: str) -> pd.DataFrame:
    """
    Adds extra kinase annotations to the dataframe.

    :param df: dataframe with measured peptide intensities
    :param extra_kinase_annot: path to file with extra kinase annotations
    """
    logger.info('Extra kinase annotation')

    extra_kinase_annot_df = pd.read_excel(extra_kinase_annot)
    df['Kinase_high_conf'] = df.index.to_series().map(extra_kinase_annot_df.set_index('Modified sequence')['PSP Kinase'])
    return df


def add_psp_urls(pp: pd.DataFrame) -> pd.DataFrame:
    """
    Function to add column to dataframe with URLs to proteins/isoforms that the
    peptides in the data belongs to:  https://www.phosphosite.org/. It uses already
    annotated information from PSP to check if any annotation exists. If it does, the URL
    is created from template and concatenated with the uniprot ID.

    :param pp: df to annotate with URL to PhosphoSitePlus
    :return: df with added annotation of URL to PhosphoSitePlus
    """
    pp[['PSP_URL', 'PSP_URL_extra']] = pp[['Start positions', 'Proteins']].apply(add_url_columns, axis=1, result_type="expand")
    return pp


def add_url_columns(row) -> Tuple[str, List]:
    start_positions, proteins = row
    # create boolean list for (p-site, protein) pairs found in PSP or not
    # check for any modified peptide with start position different from -1
    found_psites = [int(value) > 0 for value in start_positions.split(';') if value != '']
    # row[0] is integer index of row and row[1] is column value
    proteins = list(compress(proteins.split(';'),
                             found_psites))

    URLs = list()
    main_url = ""
    main_found = False

    # There can be found more than one main protein URL but then the first is used
    # and the rest goes with the isoform URLs
    url_start = "https://www.phosphosite.org/uniprotAccAction?id="
    for index, protein in enumerate(proteins):
        # don't allow isoforms (recognizable by "-" in their protein IDs) as main protein
        if '-' not in protein and not main_found:
            main_url = "=HYPERLINK(\"" + str(url_start) + str(protein) + "\")"
            main_found = True
        else:
            URLs.append(str(url_start) + str(protein))

    return main_url, URLs


def prot_clinical_annotation(df: pd.DataFrame,
                           annot_files: Union[str, Path],
                           data_type: str,
                           annot_type: str) -> Tuple[pd.DataFrame, Dict]:
    """
    Adds columns with basket annotations and weights to a dataframe

    :param df: dataframe with a 'Gene names' column to be annotated
    :param prot_baskets: list of path(s) to file(s) with annotations
    :param data_type: either 'pp' for phospho or 'fp' for full proteome
    :param basket_type: either 'basket', 'other', 'rtk' corresponding to the sheets in the annotation file
    """
    logger.info(f'Proteomics baskets annotation {data_type} {annot_type}')

    # Some dataframes might have empty cells so let's exchange them with nans
    df = df.replace(r'^\s*$', np.nan, regex=True)
    annot_dict = read_clinical_annotation(annot_files, data_type, annot_type)

    if 'fp' in data_type:
        gene_df = df.index.to_frame()
    elif 'pp' in data_type:
        gene_df = df[['Gene names']].fillna("")

    if 'POI' in annot_type:

        df[annot_type] = gene_df.apply(map_identifier_list_to_annot_types, annot_dict=annot_dict,
                                                                annot_type=annot_type,
                                                                with_weights=False,
                                                                axis=1)
    else:
        df[[annot_type, f'{annot_type}_weights']] = gene_df.apply(map_identifier_list_to_annot_types, annot_dict=annot_dict,
                                                                    annot_type=annot_type,
                                                                    with_weights=True,
                                                                    axis=1, result_type="expand")
            
    return df, annot_dict


def map_identifier_list_to_annot_types(identifier_list: pd.Series,
                                   annot_dict: Dict[str, str],
                                   annot_type: str,
                                   with_weights: bool) -> pd.Series:
    """
    Takes a list of semicolon separated identifiers and returns the annot_levels
    matching the first identifier with annotations and weights
    Input identifier list has to be pd.Series and if method used via apply it has to be of type dataframe
    """
    # TODO: throw error if no annot_dict given
    # TODO: make less hardcoded and optimize 

    annotations = []
    for identifier in identifier_list[0].split(';'):

        # TODO: instead of else if statements do small functions?
        if identifier in annot_dict:

            # should we have two dictionaries or should we split the dictionary output in case there is more than one group?
            annot_type_in_column = annot_type
            if 'POI' in annot_type:
                annot_type_in_column = 'TOPAS_subscore'

            groups = annot_dict[identifier]['GROUP'].split(';')
            annot_group = annot_dict[identifier][annot_type_in_column].split(';')
            annot_weight = annot_dict[identifier]['weight'].split(';')


            for i, group in enumerate(groups):
                # for POI only attempt dict when group is OTHER
                if group == 'OTHER' and 'POI' in annot_type:
                    annotations.append(annot_group[i])
                    
                # for TOPAS score annotations only attempt dict when group is not OTHER
                elif group != 'OTHER' and 'POI' not in annot_type:

                    if with_weights:
                        annotations.append([annot_group[i], annot_weight[i]])
                    else:
                        annotations.append(annot_group[i])
                else:
                    continue

    if 'POI' in annot_type:
        return ';'.join(annotations)
    else:
        # TODO: use a function please
        if with_weights:
            if len(annotations) > 1:
                group_names, weights = zip(*annotations)
                # take set of group names and join them with ;
                group_names = ";".join(set(group_names))
                # take set of weights and join them with ;
                weights = ";".join(weights)
                annotations = [group_names, weights]
            else:
                annotations = annotations[0] if annotations else ['']
        else:
            if len(annotations) > 0:
                annotations = ";".join(annotations)

        return pd.Series(annotations, dtype="object")
      

def create_identifier_to_basket_dict(basket_annotation: pd.DataFrame, data_type: Union[str, None] = 'fp',
                                     identifier_type: str = 'gene') -> Dict[str, str]:
    """
    collect all the baskets per gene in a dictionary of {'gene_name': 'basket1;basket2;...'}
    """
    basket_annotation = basket_annotation[basket_annotation['GROUP'] != 'OTHER']
    if 'fp' in data_type:
        accepted_type = ['expression', 'kinase activity']
        # remove non fp types
        basket_annotation = basket_annotation[basket_annotation['LEVEL'].isin(accepted_type)]
        basket_annotation = basket_annotation.groupby([identifier_type]).agg(lambda x: ";".join(map(str, x)))
    elif 'pp' in data_type:
        accepted_type = ['phosphorylation', 'important phosphorylation', 'kinase activity']
        # remove non pp types
        basket_annotation = basket_annotation[basket_annotation['LEVEL'].isin(accepted_type)]
        basket_annotation = basket_annotation.groupby([identifier_type]).agg(lambda x: ";".join(map(str, x)))
    annot_dict = basket_annotation.to_dict('index')
    return annot_dict


def create_identifier_to_poi_dict(basket_annotation: pd.DataFrame, data_type: Union[str, None] = 'fp',
                                  identifier_type: str = 'gene') -> Dict[str, str]:
    """
    """
    basket_annotation = basket_annotation[basket_annotation['GROUP'] == 'OTHER']
    basket_annotation = basket_annotation.groupby([identifier_type]).agg(lambda x: ";".join(map(str, x)))
    annot_dict = basket_annotation.to_dict('index')
    return annot_dict


def read_clinical_annotation(annot_file: str, data_type: str, annot_type: str = 'TOPAS_score') -> pd.DataFrame:
    topas_annotation = pd.read_excel(annot_file)
    topas_annotation = utils.whitespace_remover(topas_annotation)
    topas_annotation['topas_subscore_level'] = topas_annotation['TOPAS_SUBSCORE'] + " - " + topas_annotation['LEVEL']  # basket_annotation['BASKET'] + " - " +
    topas_annotation['WEIGHT'] = topas_annotation['WEIGHT'].fillna(1)  # empty cell in WEIGHT column means weight = 1
    topas_annotation = topas_annotation.rename(
        {'TOPAS_SCORE': 'TOPAS_score', 'TOPAS_SUBSCORE': 'TOPAS_subscore', 'WEIGHT': 'weight', 'GENE NAME': 'gene'}, axis=1)
    
    if 'POI' in annot_type:
        return create_identifier_to_poi_dict(topas_annotation, data_type)
    elif annot_type == 'TOPAS_score':
        return create_identifier_to_basket_dict(topas_annotation, data_type)
    else:
        return create_identifier_to_basket_dict(topas_annotation, data_type)


"""
python3 -m bin.clinical_tools -c config_patients.json -i <input_tsv> -o <output_tsv>
"""
if __name__ == "__main__":
    import argparse
    import json
    import time

    # Parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", required=True,
                        help="Absolute path to configuration file.")
    parser.add_argument("-i", "--input_file", required=True,
                        help="Absolute path to a tab separated file.")
    parser.add_argument("-o", "--output_file", required=True,
                        help="Absolute path to output file.")
    parser.add_argument("-t", "--data_type", default='fp',
                        help="Data type, either 'pp' or 'fp' (default: 'fp').")

    args = parser.parse_args()

    with open(args.config, "r") as f:
        configs = json.load(f)

    index_col = 'Gene names'
    if args.data_type == 'pp':
        index_col = 'Modified sequence'

    df = pd.read_csv(args.input_file, sep='\t', index_col=index_col)
    basket_file = configs["clinic_proc"]["prot_baskets"]

    # Start pipeline
    t0 = time.time()
    # TODO: fix.. outdated basket scheme
    for basket_type, baskets in zip(['basket', 'rtk'], [basket_file, basket_file]):
        df = prot_clinical_annotation(df, baskets,
                                    data_type=args.data_type,
                                    basket_type=basket_type)

    df.to_csv(args.output_file, sep='\t')

    t1 = time.time()
    total = t1 - t0
    logger.info(f"Basket annotation finished in {total} seconds")
