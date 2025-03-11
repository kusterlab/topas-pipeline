# %%
import os
import sys
import urllib
import argparse

import pandas as pd

import bin.config as config
import bin.TOPAS_scoring_functions as scoring

pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)


# TODO: Use pathway files from TCGA paper

def get_keggdicts():
    url_response = urllib.request.urlopen('https://rest.kegg.jp/conv/uniprot/hsa')
    uniprot_to_kegggene = url_response.read().decode(url_response.headers.get_content_charset())
    uniprot_to_kegggene = uniprot_to_kegggene.split('\n')
    for entry in uniprot_to_kegggene:
        if len(entry) <= 2:  # Filtering out entries with 2 or fewer chars; empty lines filtered here
            uniprot_to_kegggene.remove(entry)
    uniprot_to_kegggene = {splitstring.split('\t')[1].split(':')[1]: splitstring.split('\t')[0] for splitstring in
                           uniprot_to_kegggene}

    url_response = urllib.request.urlopen('https://rest.kegg.jp/link/hsa/pathway')
    kegggene_to_keggpw = url_response.read().decode(url_response.headers.get_content_charset())
    splitstring = kegggene_to_keggpw.split('\n')
    for entry in splitstring:
        if len(entry) <= 2:
            splitstring.remove(entry)
    kegggene_to_keggpw = {key: [] for key in list({i.split('\t')[1] for i in splitstring})}
    for entry in splitstring:
        entrysplit = entry.split('\t')
        kegggene_to_keggpw[entrysplit[1]].append(entrysplit[0])

    url_response = urllib.request.urlopen('https://rest.kegg.jp/list/pathway/hsa')
    keggpw_to_pathway = url_response.read().decode(url_response.headers.get_content_charset())
    keggpw_to_pathway = keggpw_to_pathway.split('\n')
    for entry in keggpw_to_pathway:
        if len(entry) <= 2:
            keggpw_to_pathway.remove(entry)
    keggpw_to_pathway = {splitstring.split('\t')[0]: splitstring.split('\t')[1] for splitstring in keggpw_to_pathway}

    del entry, entrysplit, splitstring, url_response

    return uniprot_to_kegggene, kegggene_to_keggpw, keggpw_to_pathway


def annotate_kegg_pathways(patient_df: pd.DataFrame, uniprot_to_kegggene: dict, kegggene_to_keggpw: dict,
                           keggpw_to_pathway: dict):
    """

    :param patient_df:
    :param uniprot_to_kegggene:
    :param kegggene_to_keggpw:
    :param keggpw_to_pathway:
    :return: patient_df with annotated pathway column
    """
    # 1. Add pathway column by copying protein column
    patient_df['Pathway'] = patient_df['Proteins'].copy()

    # 1.2 Explode; 1 protein per peptide
    patient_df['Pathway'] = patient_df['Pathway'].str.split(';')
    patient_df = patient_df.explode('Pathway')

    # 1.3 Filter out isoforms
    patient_df = patient_df[~patient_df['Pathway'].str.contains(r'-\d+$')]

    # 2. Replacement in pathway column; uniprot to kegg gene. Filter out NaNs
    #    pd.series.replace() is extremely slow here; using map instead
    patient_df['Pathway'] = patient_df['Pathway'].map(uniprot_to_kegggene.get)
    patient_df = patient_df[~patient_df['Pathway'].isna()]

    # 3. List builder; replace kegg gene with list of all kegg pathways the gene is associated with.
    patient_df['Pathway'] = patient_df['Pathway'].map(kegggene_to_keggpw.get)
    patient_df = patient_df[~patient_df['Pathway'].isna()]

    # 3.2 Explode; 1 kegg pathway per peptide
    patient_df = patient_df.explode('Pathway')

    # 4. Replace with human-readable names
    patient_df['Pathway name'] = patient_df['Pathway'].map(keggpw_to_pathway.get)
    patient_df['Pathway name'] = patient_df['Pathway name'].str.replace(r'- Homo sapiens \(human\)$', '', regex=True)
    return patient_df


if __name__ == '__main__':

    # comment: too long

    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", required=True, default='/home/fhamood/PycharmProjects/WP3_Pipeline/wp3_sample_pipeline/config_patients.json',
                        help="Absolute path to configuration file.")
    args = parser.parse_args(sys.argv[1:])

    configs = config.load(args.config)

    if not os.path.exists(configs["results_folder"] + '/pathway_results'):
        os.makedirs(configs["results_folder"] + '/pathway_results')

    # if not isinstance(preprocessed_df, pd.DataFrame):
    preprocessed_df = scoring.topas_score_preprocess(configs["results_folder"])

    pathway_df = scoring.calculate_psite_weights(preprocessed_df)

    uniprot_to_kegg_gene, kegg_gene_to_kegg_pathway, kegg_pathway_to_human_readable = get_keggdicts()
    pathway_df = annotate_kegg_pathways(pathway_df, uniprot_to_kegg_gene, kegg_gene_to_kegg_pathway,
                                        kegg_pathway_to_human_readable)

    pathway_summed_weights = scoring.calculate_modified_sequence_weights(pathway_df, 'Pathway name')
    pathway_capped_values = scoring.cap_zscores_and_weights(pathway_summed_weights)

    # %% calculate weighted peptide z-Scores
    pathway_scored_peptides = scoring.calculate_weighted_z_scores(pathway_capped_values)
    pathway_scored_peptides.to_csv(configs["results_folder"] + '/pathway_results/scored_peptides.tsv', sep='\t', index=False)

    # %% groupby drug and calculate pathway score
    pathway_first_level_scores = scoring.sum_weighted_z_scores(pathway_scored_peptides, by='Pathway name')

    # %% histograms for checking normal looking distribution
    scoring.plot_histograms_to_check_normality(pathway_first_level_scores)

    # %% 2nd level z-Scores
    pathway_scores = scoring.second_level_z_scoring(pathway_first_level_scores, 'Pathway name')
    pathway_scores.to_csv(configs["results_folder"] + '/pathway_results/pathway_scores.tsv', sep='\t', index=False)

