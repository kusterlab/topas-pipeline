import os
import re
import warnings
import logging
from typing import Dict, List

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib
import ast
from job_pool import JobPool

matplotlib.use('pdf')
import matplotlib.pyplot as plt

import bin.basket_scoring as basket_scoring

logger = logging.getLogger(__name__)


def create_plots_and_tables(results_folder: str, data_types: List[str]):
    # create basket score plots (~7 minutes)
    basket_swarmplots(results_folder, data_types)

    # Homologous recombination deficiency (HRD) testing (~1 minute)
    # hrd_score_swarmplots(results_folder, sample_annotation, hrd_tsg_file, imputation)

    # Homologous recombination deficiency (HRD) and tumor suppressor genes (TSG) tables (<1 minute)
    # hrd_tsg_tables(results_folder, hrd_tsg_file)

    # box plots of proteins of interest (<1 minute)
    # protein_plots(results_folder, poi_file)


def hrd_score_swarmplots(results_folder, sample_annotation, hrd_tsg_file, imputation, num_threads=4):
    logger.info('\nCalculate HRD score and create swarmplots\n')
    fp_zscore = pd.read_csv(os.path.join(results_folder, 'full_proteome_measures_z.tsv'), sep='\t', index_col='Gene names')

    sample_annotation = pd.read_csv(sample_annotation)
    sample_annotation = sample_annotation[sample_annotation['Failed'] != 'x']

    hrd_tsg = pd.read_csv(hrd_tsg_file)

    df2 = fp_zscore.copy()
    for row_index in range(0, fp_zscore.shape[0]):
        if fp_zscore.index[row_index] not in hrd_tsg['HRD'].values.tolist():
            df2 = df2.drop(fp_zscore.index[row_index])

    # replace nan with -1.5 to also weight missing TSG
    if imputation == "FALSE":
        df2 = df2.replace(np.nan, -1.5)

    # summing over patients
    summed_hrd = df2.rename(columns=dict(zip(df2.columns, sample_annotation['Sample name'])))
    summed_hrd = summed_hrd.sum(axis=0)

    z_scores = pd.DataFrame({'Sample': summed_hrd.index,
                             'Z-score': summed_hrd.values})
    z_scores = pd.melt(z_scores, id_vars='Sample',
                       value_vars=['Z-score'])

    plot_folder = os.path.join(results_folder, 'hrd_tsg')
    create_and_save_swarmplots_parallel(create_and_save_hrd_plot, plot_folder, z_scores, num_threads, summed_hrd.index)


def hrd_tsg_tables(results_folder, hrd_tsg_file):
    logger.info('Write HRD + TSG tables')
    fp_rank = pd.read_csv(os.path.join(results_folder, 'full_proteome_measures_rank.tsv'), sep='\t', index_col='Gene names')
    fp_zscore = pd.read_csv(os.path.join(results_folder, 'full_proteome_measures_z.tsv'), sep='\t', index_col='Gene names')

    hrd_tsg = pd.read_csv(hrd_tsg_file)

    path = os.path.join(results_folder, 'hrd_tsg_tables')
    os.makedirs(path, mode=0o700, exist_ok=True)

    poi_types = [hrd_tsg['HRD'], hrd_tsg['TSG']]
    for poi_type in ['HRD', 'TSG']:
        logger.info(poi_type)

        protein_of_interest_df = hrd_tsg[[poi_type]].set_index(poi_type)
        fp_rank_filtered = protein_of_interest_df.join(fp_rank, how='left')
        fp_zscore_filtered = protein_of_interest_df.join(fp_zscore, how='left')

        patients = [col_name.replace('zscore_', '') for col_name in fp_zscore_filtered.columns if col_name.startswith('zscore_')]

        for patient in patients:
            patient_df = pd.concat([fp_rank_filtered[f'rank_{patient}'], fp_zscore_filtered[f'zscore_{patient}']], axis=1)
            patient_df.columns = ['Rank', 'Z-score']
            patient_df.to_csv(os.path.join(results_folder, 'hrd_tsg_tables', f'{str(poi_type)}_fp_{patient}.csv'))


def protein_plots(results_folder, poi_file):
    fp = pd.read_csv(os.path.join(results_folder, 'annot_fp.csv'), index_col='Gene names')
    pp = pd.read_csv(os.path.join(results_folder, 'annot_pp.csv'), index_col='Gene names')
    fp_rank = pd.read_csv(os.path.join(results_folder, 'full_proteome_measures_rank.tsv'), sep='\t', index_col='Gene names')
    fp_zscore = pd.read_csv(os.path.join(results_folder, 'full_proteome_measures_z.tsv'), sep='\t', index_col='Gene names')

    fp = fp.filter(regex=r'[A-Z]\d+-')
    pp = pp.set_index(['Modified sequence', pp.index])
    pp = pp.filter(regex=r'[A-Z]\d+-')

    # read in proteins of interest
    poi = pd.read_csv(poi_file, sep='\t')
    if 'FP_group' in poi.columns:
        fp_poi = poi[['FP', 'FP_group']]
        fp_dict = dict(zip(fp_poi['FP'], fp_poi['FP_group']))
    else:
        fp_poi = poi['FP']
    if 'PP_group' in poi.columns:
        pp_poi = poi[['PP', 'PP_group']]
        pp_dict = dict(zip(zip(pp_poi['PP'], pp_poi['PP_group']), pp_poi['PP_group']))
    else:
        pp_poi = poi['PP']
    fp_poi, pp_poi = fp_poi.dropna(), pp_poi.dropna()

    path = os.path.join(results_folder, 'boxplots_POI')
    os.makedirs(path, mode=0o700, exist_ok=True)

    # protein of interest boxplot
    all_plots_fp = {}
    all_plots_fp['Patients'] = fp.columns.values.tolist()
    for row_index in range(0, fp.shape[0]):
        if fp.index[row_index] in fp_poi.values:
            all_plots_fp[str(fp.index[row_index])] = fp.iloc[row_index, :].values.tolist()

    all_plots_pp = {}
    all_plots_pp['Patients'] = pp.columns.values.tolist()
    if 'PP_group' in poi.columns:
        for row_index in range(0, pp.shape[0]):
            if pp.index[row_index] in pp_dict.keys():
                all_plots_pp[str(pp.index[row_index])] = pp.iloc[row_index, :].values.tolist()

    boxplot_intensity_fp, boxplot_intensity_pp = pd.DataFrame(all_plots_fp), pd.DataFrame(all_plots_pp)
    boxplot_intensity_fp = pd.melt(boxplot_intensity_fp, id_vars='Patients')
    boxplot_intensity_pp = pd.melt(boxplot_intensity_pp, id_vars='Patients')

    def f(x):
        try:
            llist = ast.literal_eval(x)
            llist = [string.strip() for string in llist]
            return llist
        except Exception:
            return x

    if 'FP_group' in poi.columns:
        for index in boxplot_intensity_fp.index:
            boxplot_intensity_fp.loc[index, 'Group'] = fp_dict[boxplot_intensity_fp.loc[index, 'variable']]

    # if 'PP_group' in poi.columns:
    #     for index in boxplot_intensity_pp.index:
    #         boxplot_intensity_pp.loc[index, 'Group'] = pp_dict[boxplot_intensity_pp.loc[index, 'variable']]

    boxplot_intensity_pp['variable'] = boxplot_intensity_pp['variable'].apply(f)
    boxplot_intensity_pp['Group'] = boxplot_intensity_pp['variable'].apply(lambda x: x[1])
    boxplot_intensity_pp['variable'] = boxplot_intensity_pp['variable'].apply(lambda x: x[0])

    # check if patient specific:
    patient_match = re.search('[A-Z]\d+-[A-Z|\d]+-[A-Z|\d]+', poi_file)
    if patient_match:
        patient_of_interest = patient_match.group(0)

    # if entity specific
    # entity = True
    # RMS entity
    # # Todo make as input
    # # entity = ['I076-014-129774', 'I054-033-226658', 'I003-006-103006', 'I024-015-104566', 'I024-034-189638', 'I028-003-84352',
    # #           'I034-044-186620-R2', 'I034-052-222996', 'I062-008-87120-R2', 'I078-024-223164', 'I010-021-226690',
    # #           'I002-010-106444-R2',
    # #           'I007-007-80690-R2', 'I007-020-1007541', 'I007-031-108742', 'I007-039-130734', 'I015-006-107208', 'I022-018-98376',
    # #           'I024-020-127332-R2', 'I036-012-107936-R2', 'I043-001-80842', 'I043-005-130270', 'I043-005-95540',
    # #           'I052-003-93894-R2',
    # #           'I137-003-84050', 'H021-GDTFYK-M2-R2', 'H021-LGAPJC-M1', 'I010-022-226744']
    #
    # if entity:
    #     boxplot_intensity_fp = boxplot_intensity_fp[boxplot_intensity_fp['Patients'].isin(entity)]
    #

    for i, group in enumerate(boxplot_intensity_pp.groupby('Group')):
        plt.figure(figsize=(10, 12))
        group_pp = group[1]

        for j, patient in enumerate(all_plots_pp['Patients']):  # loop not necessary
            if patient == patient_of_interest:  # add 'zscore_'+ in front if using that data
                temp_patient = group_pp[group_pp['Patients'] == patient]

                ax = sns.boxplot(x="variable", y="value", data=group_pp, color='lightgrey',
                                 order=pp_poi[pp_poi['PP_group'] == group[0]]['PP'].unique())
                ax = sns.swarmplot(x="variable", y="value", data=group_pp,
                                   color=".25", order=pp_poi[pp_poi['PP_group'] == group[0]]['PP'].unique())

                ax = sns.swarmplot(x="variable", y="value", data=temp_patient,
                                   color="red", size=10, order=pp_poi[pp_poi['PP_group'] == group[0]]['PP'].unique())

                ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
                ax.set_ylabel("Intensity", fontsize=16)
                ax.set(xlabel=None)
                title = patient + ', PP'
                plt.title(title, fontsize=18)
                locs, labels = plt.xticks()
                plt.setp(labels, rotation=45)

                plt.tight_layout()
                plt.savefig(os.path.join(results_folder, 'boxplots_POI', f'{patient}_PP_{group[0]}_POI_intensity.png'), dpi=300)
                plt.close()

    # add FP or PP
    # for when there is a group
    for i, group in enumerate(boxplot_intensity_fp.groupby('Group')):
        plt.figure(figsize=(10, 12))
        group_fp = group[1]

        for j, patient in enumerate(all_plots_fp['Patients']):  # loop not necessary
            if patient == patient_of_interest:  # add 'zscore_'+ in front if using that data
                temp_patient = group_fp[group_fp['Patients'] == patient]
                logger.info(fp_poi[fp_poi['FP_group'] == group[0]]['FP'].unique())
                ax = sns.boxplot(x="variable", y="value", data=group_fp, color='lightgrey',
                                 order=fp_poi[fp_poi['FP_group'] == group[0]]['FP'].unique())
                ax = sns.swarmplot(x="variable", y="value", data=group_fp,
                                   color=".25", order=fp_poi[fp_poi['FP_group'] == group[0]]['FP'].unique())

                ax = sns.swarmplot(x="variable", y="value", data=temp_patient,
                                   color="red", size=10, order=fp_poi[fp_poi['FP_group'] == group[0]]['FP'].unique())

                ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
                ax.set_ylabel("Intensity", fontsize=16)
                ax.set(xlabel=None)
                title = patient + ', FP'
                plt.title(title, fontsize=18)
                locs, labels = plt.xticks()
                plt.setp(labels, rotation=45)

                plt.tight_layout()
                plt.savefig(os.path.join(results_folder, 'boxplots_POI', f'{patient}_{group[0]}_POI_intensity.png'), dpi=300)
                plt.close()


def basket_swarmplots(results_folder, data_types, num_threads=4):
    """

    :param results_folder:
    """

    basket_scores = basket_scoring.read_basket_scores(results_folder, main_basket_only=False)
    all_samples = basket_scores.index

    # TODO: test old code again (is at least very close to working or perhaps even working)
    # TODO: and make this work on both generations
    if basket_scores.columns.str.contains('Scores').any():

        clinical_baskets, immune_rtk_baskets, rtk_baskets = get_baskets_order(basket_scores.columns)

        rtk_number_table = create_num_rtk_annotation_df(basket_scores, data_types)
        basket_types = [('baskets', 'Scores', clinical_baskets, None),
                        ('immune_rtk', 'Scores', immune_rtk_baskets, None),
                        ('rtk_baskets', 'RTK Scores', rtk_baskets, rtk_number_table)]

        for data_type in data_types:
            if data_type == 'pp':
                if 'Immunotherapy' in clinical_baskets:
                    clinical_baskets.remove('Immunotherapy')
                    clinical_baskets.remove('DRD')

            for basket_type, sheet_suffix, baskets, number_table in basket_types:
                if len(baskets) == 0:
                    logger.info(f"Skipping {basket_type} {data_type} box plots as no baskets were given")
                    continue

                logger.info(f"Plotting {basket_type} {data_type} box plots")
                plot_folder = os.path.join(results_folder, basket_type)

                df = create_boxplot_df(basket_scores, f'{data_type.upper()} - {sheet_suffix}', baskets)
                create_and_save_swarmplots_parallel(create_and_save_basket_boxplot, plot_folder, df, num_threads, all_samples, data_type, number_table, basket_type)
                logger.info("Plotting protein/p-site numbers box plots")
                number_table = create_num_identifications_and_annotations_df(basket_scores, data_types)
    else:
        logger.info("Plotting protein/p-site numbers box plots")
        number_table = save_num_identifications_and_annotations_df_4th_gen(basket_scores, data_types)
    number_table.to_csv(os.path.join(results_folder, 'ident_annot_numbers_new.tsv'))

    # plot_folder = os.path.join(results_folder, 'ident_annot')
    # create_and_save_swarmplots_parallel(create_and_save_id_number_boxplot, plot_folder, number_table, num_threads, all_samples, data_types)


def create_and_save_swarmplots_parallel(plot_function, plot_folder, df, num_threads, all_samples, *args):
    os.makedirs(plot_folder, mode=0o700, exist_ok=True)

    ax = create_boxplot(df)

    processingPool = JobPool(processes=num_threads)
    for patient in all_samples:
        processingPool.applyAsync(plot_function, [plot_folder, df, ax, patient, *args])

    processingPool.checkPool(printProgressEvery=num_threads)


def create_boxplot_df(basket_scores, sheet_name, baskets):
    basket_dict = {sheet_name + '.' + k: k for k in baskets}
    return create_melted_df(basket_scores, basket_dict)


def create_num_rtk_annotation_df(basket_scores, data_types):
    basket_dict = {f'{data_type.upper()} - RTK Scores.num_annotated': f'Number proteins {data_type} rtk' for data_type in data_types}
    return create_melted_df(basket_scores, basket_dict)


def create_num_identifications_and_annotations_df(basket_scores, data_types):
    basket_dict = {f'{data_type.upper()} - Scores.num_identified': f'Number proteins ident {data_type}' for data_type in data_types} | \
                  {f'{data_type.upper()} - Scores.num_annotated': f'Number proteins annot {data_type}' for data_type in data_types}
    return create_melted_df(basket_scores, basket_dict)


def save_num_identifications_and_annotations_df_4th_gen(basket_scores, data_types):
    basket_scores = basket_scores.loc[:, (basket_scores.columns.str.contains('num_identified')) | (basket_scores.columns.str.contains(
        'num_annotated'))]
    basket_dict = {f'{data_type}.num_identified': f'Number proteins ident {data_type}' for data_type in data_types} | \
                  {f'{data_type}.num_annotated': f'Number proteins annot {data_type}' for data_type in data_types}
    return create_melted_df(basket_scores, basket_dict)


def create_melted_df(basket_scores, basket_dict: Dict[str, str]):
    available_baskets = {k: v for k, v in basket_dict.items() if k in basket_scores.columns}
    boxplot_df = basket_scores[available_baskets.keys()]
    boxplot_df = boxplot_df.rename(columns=available_baskets).reset_index()
    boxplot_df = pd.melt(boxplot_df, id_vars='Sample', value_vars=available_baskets.values())
    return boxplot_df


def create_and_save_basket_boxplot(plot_folder, df, ax, patient, data_type, number_table, suffix):
    ax = add_patient_dot(df, ax, patient)

    title_suffix = ""
    if number_table is not None:
        temp_number = number_table[number_table['Sample'] == patient]
        protein_number = temp_number[temp_number['variable'] == f'Number proteins {data_type} rtk'][
            'value'].values
        title_suffix = "    " + 'n=' + str(protein_number[0])

    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    ax.set_ylabel("Basket score", fontsize=16)
    ax.set(xlabel=None)

    title = patient + ', ' + data_type.upper() + title_suffix
    plt.title(title, fontsize=18)
    locs, labels = plt.xticks()
    plt.setp(labels, rotation=45)

    plt.tight_layout()
    with warnings.catch_warnings():
        warnings.filterwarnings(action="ignore", category=UserWarning, message=r".*of the points cannot be placed.*")
        plt.savefig(os.path.join(plot_folder, f'{patient}_{data_type.upper()}_score_{suffix}.pdf'))
    plt.close()


def create_and_save_id_number_boxplot(plot_folder, number_table, ax, patient, data_types):
    ax = add_patient_dot(number_table, ax, patient)

    ax.set_xticklabels(
        [f'{data_type.upper()} ident.' for data_type in data_types] + [f'{data_type.upper()} annot.' for data_type in data_types])

    ax.set_ylabel("Number proteins/p-site", fontsize=12)
    ax.set(xlabel=None)

    plt.title(patient + ' - Number of identified and annotated '
                        'proteins/phospho peptides')

    plt.tight_layout()
    with warnings.catch_warnings():
        warnings.filterwarnings(action="ignore", category=UserWarning, message=r".*of the points cannot be placed.*")
        plt.savefig(os.path.join(plot_folder, f'{patient}_ident_annot_numbers.pdf'))
    plt.close()


def create_and_save_hrd_plot(plot_folder, z_scores, ax, patient):
    ax = add_patient_dot(z_scores, ax, patient)

    ax.set_ylabel("Summed z-score", fontsize=12)
    ax.set(xlabel=None)
    plt.title('HRD genes summed Z-score - ' + str(patient), fontsize=16)
    locs, labels = plt.xticks()
    plt.setp(labels, rotation=45)

    plt.tight_layout()
    with warnings.catch_warnings():
        warnings.filterwarnings(action="ignore", category=UserWarning, message=r".*of the points cannot be placed.*")
        plt.savefig(os.path.join(plot_folder, f'{patient}_summed_zscore_hrd_genes.pdf'))
    plt.close()


def create_boxplot(df):
    plt.figure(figsize=(12, 8))
    with warnings.catch_warnings():
        warnings.filterwarnings(action="ignore", category=UserWarning, message=r".*of the points cannot be placed.*")
        ax = sns.boxplot(x="variable", y="value", data=df, color='lightgrey')
        ax = sns.swarmplot(x="variable", y="value", data=df,
                           color=".25", size=4)
    return ax


def add_patient_dot(df, ax, patient):
    temp_patient = df[df['Sample'] == patient]
    ax = sns.swarmplot(x="variable", y="value", data=temp_patient,
                       color="red", size=10, ax=ax)
    return ax


def get_baskets_order(columns):
    # TODO: move this to a json file
    if 'FP - Scores.Signal Transduction' in columns:  # 2nd generation baskets
        clinical_baskets = ['NRTK',
                            'RTK',
                            'Signal Transduction',
                            'MAPK Ras',
                            'PI3K-AKT-mTOR',
                            'Immunotherapy',
                            'JNK/p38 MAPK',
                            'NFkb',
                            'ECM',
                            'Cell Cycle']
        immune_rtk_baskets = []
        rtk_baskets = ['EGFR',
                       'ERBB',
                       'PDGFRB',
                       'PDGFRA',
                       'VEGFR',
                       'FGFR',
                       'MET',
                       'ALK',
                       'KIT',
                       'NTRK1',
                       'NTRK2',
                       'NTRK3',
                       'MERTK',
                       'AXL',
                       'DDR1/2',
                       'IGFR1',
                       'EPH',
                       'ROR1/2',
                       'OSMR',
                       'AURKA/B'
                       ]

    elif 'FP - Scores.TK' in columns:  # 3rd generation baskets
        clinical_baskets = ['TK',
                            'RTK signal transmission',
                            'MAPK Ras',
                            'PI3K-AKT-mTOR',
                            'Immunotherapy',
                            'JNK/p38 MAPK',
                            'NFkb',
                            'Developmental',
                            'DRD',
                            'Cell Cycle']
        immune_rtk_baskets = []
        rtk_baskets = ['EGFR',
                       'ERBB',
                       'PDGFRB',
                       'PDGFRA',
                       'VEGFR',
                       'FGFR',
                       'MET',
                       'ALK',
                       'KIT',
                       'NTRK1',
                       'NTRK2',
                       'NTRK3',
                       'MERTK',
                       'AXL',
                       'DDR1/2',
                       'IGFR1',
                       'EPH',
                       'ROR1/2',
                       'OSMR',
                       'AURKA/B'
                       ]
    else:  # 1st generation
        clinical_baskets = ['DNA damage',
                            'MAPK RAS signaling',
                            'PI3K AKT mTOR',
                            'Cell cycle',
                            'Wnt Notch Hedgehog']
        immune_rtk_baskets = ['Immunotherapy',
                              'RTKs']
        rtk_baskets = ['ErbB',
                       'EGF/EGFR',
                       'ERBB2',
                       'ERBB4',
                       'FGFR1',
                       'FGFR2',
                       'FGFR3',
                       'MET',
                       'MST1',
                       'PDGF',
                       'PDGFRb',
                       'KIT',
                       'NTRK1',
                       'NTRK2',
                       'NTRK3',
                       'VEGF',
                       'FLT3',
                       'JAK/STAT',
                       'RET',
                       'EPH-Ephrin',
                       'IGF1R',
                       'HGF',
                       'AXL',
                       'IL6R',
                       'Oncostatin M'
                       ]
    return clinical_baskets, immune_rtk_baskets, rtk_baskets


if __name__ == '__main__':
    import argparse

    from . import config

    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", required=True,
                        help="Absolute path to configuration file.")
    args = parser.parse_args(sys.argv[1:])

    configs = config.load(args.config)

    create_plots_and_tables(configs["results_folder"], data_types=configs["data_types"])
