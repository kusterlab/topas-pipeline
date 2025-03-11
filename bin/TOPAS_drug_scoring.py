from pathlib import Path
import re
import sys
import logging
import argparse
import os
import shutil

import pandas as pd
import numpy as np
import toml

import psite_annotation as pa
import bin.config as config
import bin.TOPAS_scoring_functions as scoring

# hacky way to get the package logger instead of just __main__ when running as a module
logger = logging.getLogger(__package__ + "." + Path(__file__).stem)

pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)

mode = 3
mode_dict = {
    1: {'name': 'A431', 'cells': ['A431']},
    2: {'name': 'all_but_A431', 'cells': ['SKES1', 'A204', 'MESSA', 'SKLMS1']},
    3: {'name': 'all', 'cells': ['SKES1', 'A204', 'MESSA', 'SKLMS1', 'A431']}
}


def perform_grouped_aggregations(grouped_df):
    aggregation_functions = {
        'Site positions': pd.NamedAgg(column='Site positions', aggfunc=csv_list_unique),
        'Proteins': pd.NamedAgg(column='Proteins', aggfunc=csv_list_unique),
        'pEC50': pd.NamedAgg(column='pEC50', aggfunc=csv_list_unique),
        'Curve Regulation': pd.NamedAgg(column='Curve Regulation', aggfunc=csv_list_unique),
        'Start positions': pd.NamedAgg(column='Start positions', aggfunc=csv_list_unique),
        'End positions': pd.NamedAgg(column='End positions', aggfunc=csv_list_unique),
        'Site sequence context': pd.NamedAgg(column='Site sequence context', aggfunc=csv_list_unique),
        'Datasets': pd.NamedAgg(column='Dataset', aggfunc=csv_list_unique),
        'Dataset Count': pd.NamedAgg(column='Dataset', aggfunc='count')
    }
    if 'Dataset Count' in grouped_df.first().columns:
        aggregation_functions.pop('Datasets')
        aggregation_functions['Dataset Count'] = pd.NamedAgg(column='Dataset Count', aggfunc=max)
        aggregation_functions['Modified sequence'] = pd.NamedAgg(column='Modified sequence', aggfunc=csv_list_unique)

    grouped_df = grouped_df.agg(**aggregation_functions)
    return grouped_df


def merge_with_drug_annotations(pp_df, drug_annotations):
    pp_df = pd.merge(left=pp_df, right=drug_annotations, on='Site positions', how='inner')
    return pp_df


def filter_for_drug_psites(pp_df, drug_path, cell_lines="all"):
    drug_annotations = get_drug_annotation_dict(cell_lines, drug_path)
    pp_df = pp_df[pp_df['Site positions'].isin(drug_annotations['Site positions'])]
    # TODO: This filtering step for non-isoform peptides is not needed; performed in score_preprocess(). Remove!
    pattern = re.compile(r'-\d+_[STY]')
    non_iso_sites = [site for site in pp_df['Site positions'] if not re.search(pattern, site)]
    pp_df = pp_df.loc[pp_df['Site positions'].isin(non_iso_sites)]
    return pp_df, drug_annotations


def map_psite_to_drug(psite_df, annotations):
    for psite in psite_df[0].split(';'):
        if psite in annotations:
            return pd.Series(annotations[psite].get_drug_annotation())

    return pd.Series({'basket': "", 'weights': "", 'pEC50': "", 'regulation': ""})


def get_drug_annotation_dict(cell_lines: str, drug_path):
    """
    :param cell_lines: comma-separated string of all cell lines used, or alternatively 'all' to use all available cell
    lines
    :return: multi-level dictionary mapping pSites to their pSite and their weight for score calculation
             {<pSite>: {'basket': <Drug>, 'weight': <weight>}, <pSite>: {...}, ...}
    """
    curves_df = load_topas_wp1_data(cell_lines, drug_path)
    curves_df = set_curve_weights(curves_df)

    # ['Modified sequence', 'Proteins', 'pEC50', 'Curve Regulation', 'Drug']

    # curves_df = curves_df[['basket', 'Site positions', 'weight']]

    curves_df = curves_df[curves_df['Site positions'].replace("", np.nan).notna()]
    basket_annotation = curves_df.drop(columns=['Proteins', 'Start positions', 'End positions'])
    basket_annotation = basket_annotation.rename(columns={'Modified sequence': 'Modified sequences in curves'})
    return basket_annotation


def set_curve_weights(curves_df):
    base_weight = 0.7
    cell_line_count_weight = 0.3

    # TODO: decide what to do with upgoing curves: -1, 0 or 1
    curves_df['weight'] = curves_df['Curve Regulation'].map({'down': 1, 'up': 0})
    curves_df['weight'] *= base_weight + (cell_line_count_weight * curves_df['Dataset Count'])

    return curves_df


def map_drugs(dataset):
    # TODO: add to config w lambda function
    tomlpath = f'/media/kusterlab/internal_projects/active/TOPAS/WP11/Data_Sets/TOPAS_{dataset}/0__tomls/'
    drugmap = dict()
    for file in os.listdir(tomlpath):
        drug_toml = toml.load(tomlpath + file)
        drugmap[file.split('.')[0]] = drug_toml['Sample']['drug']
    return drugmap


def csv_list_unique(x, separator=';'):
    return separator.join(map(str, list(set(x))))


def list_unique(x):
    return list(set(x))


def aggregate_curve_regulation(x):
    # TODO: Replace with classifier of "up", "down", or "ambiguous" and change check in set_curve_weights()
    return ",".join(map(str, list(dict.fromkeys(x))))


def change_phospho_notation(mod_seq):
    return re.sub(r"([STY])\(ph\)", "p\1", mod_seq)


def remove_ambiguous_curves(curves_df):
    curves_df = curves_df.loc[curves_df['Curve Regulation'].isin(['up', 'down'])]
    return curves_df


def load_topas_wp1_data(cell_dataset, drug_path):
    # TODO: tidy up to not use mode_dict and instead use input
    # excluded_drugs = ['Anlotinib', 'Barasertib', 'Brivanib alaninate', 'Buparlisib',
    #                   'Capmatinib', 'Cediranib', 'Crenolanib', 'Infigratinib']
    excluded_drugs = []

    curves_df = pd.DataFrame()
    for cell_line in mode_dict[mode]["cells"]:
        cell_line_curves = read_regulated_curves_file(cell_line, excluded_drugs, drug_path)
        cell_line_curves['Dataset'] = cell_line
        curves_df = pd.concat([curves_df, cell_line_curves])
    curves_df = filter_bad_curves(curves_df)

    curves_df = aggregate_curves_multiple_cell_lines(curves_df)
    curves_df = remove_ambiguous_curves(curves_df)

    curves_df['Site positions'] = curves_df['Site positions'].apply(lambda x: x.split(";"))
    curves_df = curves_df.explode('Site positions')

    curves_df = curves_df.groupby(['Drug', 'Site positions'], as_index=False)
    curves_df = perform_grouped_aggregations(curves_df)
    return curves_df


def filter_bad_curves(curves_df):
    return curves_df[(curves_df['pEC50'] < 10.5) & (curves_df['pEC50'] > 5)]


def read_regulated_curves_file(cell_line, excluded_drugs, drug_path):
    cols = ['Modified sequence', 'Proteins', 'pEC50', 'Curve Regulation', 'Experiment']
    # TODO: move out of function + give base location in config?
    curvefile = f'/media/kusterlab/internal_projects/active/TOPAS/WP11/Data_Sets/TOPAS_{cell_line}/10__analysis/all_regulated_curves.txt'
    curvefile_save = os.path.join(drug_path, f'all_regulated_curves_{cell_line}.txt')
    curves_df = pd.read_csv(curvefile, usecols=cols, sep='\t')
    shutil.copyfile(curvefile, curvefile_save)

    curves_df = curves_df.replace(map_drugs(cell_line))
    curves_df = curves_df.rename(columns={'Experiment': 'Drug'})
    curves_df = curves_df[~curves_df['Drug'].isin(excluded_drugs)]
    curves_df = phospho_annot(curves_df)
    # TODO: This could filter out curves if the same drug-site-combination shows an upgoing and downgoing curve...
    curves_df = curves_df.drop_duplicates(subset=['Drug', 'Site positions'], keep='first')
    return curves_df


def phospho_annot(df) -> pd.DataFrame:
    """
    Phospho-site annotation of experimental data using in-house developed tool (MT) based mainly on Phosphosite Plus but
     also in vitro experiments from the lab of Ishihama.

    :param df: dataframe with measured peptide intensities
    :return:

    """
    # TODO: can be found in config
    pspFastaFile = "/media/kusterlab/line_functions/bioinformatics/Databases/psite_annotation_mapping/PhosphoSitePlus/Phosphosite_seq.fasta"
    df = pa.addPeptideAndPsitePositions(df, pspFastaFile, pspInput=True)

    return df


def aggregate_curves_multiple_cell_lines(curves_df):
    curgrp = curves_df.groupby(['Drug', 'Modified sequence'], as_index=False)
    curgrp = perform_grouped_aggregations(curgrp)
    return curgrp


def drug_scoring(results_folder, preprocessed_df):
    logger.info('Running drug scoring module')

    if os.path.exists(os.path.join(results_folder, 'drug_results_all', 'drug_scores.tsv')):
        logger.info(f'Drug scoring skipped - found files already processed')
        return

    drug_path = os.path.join(results_folder, f'drug_results_{mode_dict[mode]["name"]}')

    if not os.path.exists(drug_path):
        os.makedirs(drug_path)

    preprocessed_df = scoring.read_preprocessed_df(os.path.join(results_folder, 'topas_score_preprocessed.tsv'))
    drug_df, annotations = filter_for_drug_psites(preprocessed_df, drug_path, cell_lines=mode_dict[mode]["name"])
    logger.info('  Calculate p-site weights')
    drug_df = scoring.calculate_psite_weights(drug_df)

    drug_annotated_patients = merge_with_drug_annotations(drug_df, annotations)
    drug_annotated_patients.to_csv(os.path.join(results_folder, 'drug_annotated_patients_with_weights.tsv'), sep='\t',
                                   index=False)

    logger.info('  Calculate modified sequence weights')
    drug_summed_weights = scoring.calculate_modified_sequence_weights(drug_annotated_patients, 'Drug')

    # TODO: Export here; table with uncapped weights and uncapped zscores
    drug_capped_values = scoring.cap_zscores_and_weights(drug_summed_weights)

    logger.info('  Calculate weighted z-scores')
    drug_scored_peptides = scoring.calculate_weighted_z_scores(drug_capped_values)
    drug_scored_peptides.to_csv(os.path.join(drug_path, 'scored_peptides.tsv'), sep='\t', index=False)

    logger.info('  Calculate drug scores')
    drug_first_level_scores = scoring.sum_weighted_z_scores(drug_scored_peptides, by='Drug')

    # scoring.plot_histograms_to_check_normality(drug_first_level_scores)

    logger.info('  2nd level z-scoring, adding target space and writing results')
    drug_scores = scoring.second_level_z_scoring(drug_first_level_scores, 'Drug')
    # TODO: what is target space  --> no function documentation and name could be anything
    drug_spaces = scoring.get_target_space(annotated_peptides_df=drug_annotated_patients, scored_peptides_df=drug_scored_peptides,
                                           grouping_by='Drug')
    drug_scores = pd.merge(left=drug_spaces, right=drug_scores, on='Drug', how='left').sort_values(by='Drug')
    drug_scores.to_csv(os.path.join(drug_path, 'drug_scores.tsv'), sep='\t', index=False)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", required=False,
                        default='/home/fhamood/PycharmProjects/WP3_Pipeline/wp3_sample_pipeline/config_patients.json',
                        help="Absolute path to configuration file.")
    args = parser.parse_args(sys.argv[1:])
    configs = config.load(args.config)

    preprocessed_path = os.path.join(configs["results_folder"], 'topas_score_preprocessed.tsv')
    if os.path.exists(preprocessed_path):
        print('Found preprocessed file')
        preprocessed_df = pd.read_csv(preprocessed_path, sep='\t')
    else:
        preprocessed_df = scoring.topas_score_preprocess(configs["results_folder"])

    drug_scoring(configs["results_folder"], preprocessed_df)
