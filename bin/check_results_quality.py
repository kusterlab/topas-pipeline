import os
import re
import sys
import argparse
import logging
from pathlib import Path
from collections import defaultdict
from typing import List, Dict, Tuple, Any, Union

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
from scipy import stats

import bin.config as config
from bin.qc.peptide_heatmap import plot_peptide_heatmaps
import bin.preprocess_tools as prep
import bin.pca_plotting as pca_tools
from bin import sample_metadata, sample_annotation

if sys.platform.startswith('linux'):
    matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages

logger = logging.getLogger(__name__)

plt.rcParams['text.usetex'] = False


def check_results(results_folder: Union[str, Path],
                  new_batch: int,
                  sample_annotation_file: str,
                  meta_data: str,
                  data_types: List[str]) -> None:
    # TODO: fix that this also works for only full proteome
    # TODO: make it class based?

    if len(data_types) < 2:
        return
    # create subdirectory for QC results
    os.makedirs(os.path.join(results_folder, 'Results_investigation_plots'), mode=0o700, exist_ok=True)

    # Get sample and metadata annotation
    meta_data_df = sample_metadata.load(meta_data)
    meta_data_df = pca_tools.whitespace_remover(meta_data_df)
    # sample_annotation_df = pd.read_csv(os.path.join(results_folder, 'sample_annot_filtered.csv'))
    sample_annotation_df = pd.read_csv(sample_annotation_file)
    sample_annotation_df = sample_annotation.filter_sample_annotation(sample_annotation_df, remove_qc_failed=True, remove_replicates=True)
    batches = sample_annotation_df['Batch Name'].unique()

    # Get pre- and fully processed data
    data_dict = read_patient_data(results_folder, sample_annotation_df, data_types)

    data_types = ['pp']

    # PCA    # can one loop through pct and save plots?
    plot_types = ['Intensity']
    selected_proteins = []
    principal_dict = {}
    principal_var_dict = {}
    # for data_type in data_types:
    #     principal_dfs, principal_variances, meta_data_df = pca_tools.do_pca(selected_proteins, results_folder, plot_types,
    #                                                                         sample_annotation_file,
    #                                                                         meta_data, data_type,
    #                                                                         include_replicates=True, plot=False)
    #     principal_dict[data_type] = principal_dfs
    #     principal_var_dict[data_type] = principal_variances
    #
    # # Median vs count
    # median_vs_count_dfs = get_z_score_idents([data_dict['full_proteome_measures_z'], data_dict['phospho_measures_z']],
    #                                          meta_data_df)

    # QC of ref channels
    for data_type in data_types:
        ref_channel_df = ref_channel_analysis(data_dict, batches, data_type, results_folder)

    #TODO: check that it still works and give info as return to use in plot qc report

    # plot_qc_results(results_folder, median_vs_count_dfs, principal_dict, principal_var_dict, new_batch, ref_channel_df)

    # QC report
    # create_qc_report(data_dict, meta_data_dicts)

    # Plot PCA with metadata labels
    # metadata_columns = ['Replicates', 'Batch_No', 'Sarcoma Subtype', 'Program', 'Tumor cell content']
    # meta_data_dicts = {}
    # for metadata in metadata_columns:
    #     meta_data_dicts[metadata] = get_sample_dict(sample_annotation_df, meta_data_df, remove_qc_failed=True, data_type=metadata)
    #     plot_metadata_pcas(results_folder, meta_data_dicts[metadata], metadata, data_dict)

    # quick sub cohorts plot - misdiagnosis
    # Subtype pca - chordoma (not cohort but type metadata annot)
    # add label of outlier samples.. give input
    # fp_df = data_dict['preprocessed_fp']
    # pp_df = data_dict['preprocessed_pp']
    # fp_df = fp_df.loc[:, fp_df.columns.isin(meta_data_dicts['Sarcoma Subtype'].keys())]
    # pp_df = pp_df.loc[:, pp_df.columns.isin(meta_data_dicts['Sarcoma Subtype'].keys())]
    # labels = ['I005-017-1T1-P1', 'H021-HQC7AD-T1']
    # data_dict['preprocessed_fp_ref'] = data_dict['preprocessed_fp_ref'].loc[:, ~data_dict['preprocessed_fp_ref'].columns.isin(
    #     data_dict['preprocessed_fp_ref'].filter(regex='metadata').columns.tolist())]
    # pca_colors, pca_markers = get_pca_color_marker(fp_df, data_dict['preprocessed_fp_ref'])
    # pca_colors = ['green', 'red']
    # pca_objects_list, _ = get_pca_objects([fp_df, pp_df], pct_threshold=[0.5],
    #                                       do_pca=False)
    # targets = get_pca_targets(fp_df, data_dict['preprocessed_fp_ref'])

    # principal_dfs = get_pca_dfs([fp_df, pp_df], pca_objects_list, targets, pca_type='ppca')
    # for df in principal_dfs:
    #     # replace anything but rhabdomyosarcoma
    #     df['Labels'] = df['Sample']
    #     for i, element in enumerate(df['Metadata label']):
    #         if element not in ['Rhabdomyosarcoma', 'Osteosarcoma']:
    #             df.loc[df.index[i], 'Metadata label'] = 'Other'
    #     for i, name in enumerate(df['Sample']):
    #         if name not in labels:
    #             df.loc[df.index[i], 'Labels'] = ' '
    # pc_variances = []
    # for data_type in pca_objects_list:
    #     pc_variances.append(data_type.var_exp)
    # fig, ax = plt.subplots(2, 1, figsize=(8.27, 11.69))
    # fig = plot_pca(fig, ax, df=[principal_dfs[0], principal_dfs[1]],
    #                exp_var=pc_variances, target='Metadata label', subplot=True,
    #                colors=pca_colors, markers=pca_markers, labels='Labels')  # [pc_variances[0], pc_variances[1]]
    # plt.tight_layout()
    # plt.savefig(f'{results_folder}/Results_investigation_plots/pca_patient_misdiagnosis.pdf')

    # TODO: Add basket score plots... with and without ref channels

    # Plot PCA with reference channels
    # if 'preprocessed_fp_with_ref' in data_dict.keys():
    #     plot_ref_pca(results_folder, meta_data_df, sample_annotation_df, data_dict)
    #
    #     # peptide heatmap for CDK2NA
    #     if debug:
    #         plot_peptide_heatmaps(results_folder, sample_annotation_file)

    # TODO: add to QC report
    # Preprocessing PCAs
    # preprocess_data = [data_dict['debug_preprocessed_fp2_complete_raw'],
    #                    data_dict['debug_preprocessed_fp2_before_ms1_correction'],
    #                    data_dict['debug_preprocessed_fp2_after_ms1_correction']]
    # for i, df in enumerate(preprocess_data):
    #     df = df.loc[:, df.columns.isin(meta_data_dicts['Batch_No'].keys())]
    #     pca_colors, pca_markers = get_pca_color_marker(df, meta_data_dicts['Batch_No'])
    #     pca_objects_list, _ = get_pca_objects([df], pct_threshold=[0.5], do_pca=False)
    #     targets = get_pca_targets(df, meta_data_dicts['Batch_No'])
    #     principal_dfs = get_pca_dfs([df], pca_objects_list, targets, pca_type='ppca')
    #     pc_variances = []
    #     for data_type in pca_objects_list:
    #         pc_variances.append(data_type.var_exp)
    #     fig, ax = plt.subplots(1, 1, figsize=(12, 12))
    #     fig = plot_pca(fig, ax, df=principal_dfs,
    #                    exp_var=pc_variances, target='Metadata label', subplot=False,
    #                    colors=pca_colors, markers=pca_markers)
    #     plt.tight_layout()
    #     plt.savefig(os.path.join(results_folder, 'Results_investigation_plots', f'pca_prep_{i}_pp.pdf'))


def ref_channel_analysis(data_dict, batches, data_type, results_folder):
    in_batch_overlaps = {}
    in_batch_counts = {}
    between_batches_overlap = {}
    ref_channels = get_ref_channels(data_dict[f'annot_{data_type}_with_ref'])
    intensity_channels = get_intensity_channels(data_dict[f'annot_{data_type}_with_ref'])
    intensity_range = [intensity_channels.min(axis=1), intensity_channels.max(axis=1)]
    median_intensity = intensity_channels.median(axis=1)
    std_intensity = intensity_channels.std(axis=1)
    del intensity_channels
    annot_channels = get_annot(data_dict[f'annot_{data_type}_with_ref'])
    # use regex to get batches and filter out nan to get overlap
    for batch in batches:
        temp_ref_channels = ref_channels.filter(regex=fr'Batch{batch}$')
        temp_ref_channels_filtered = temp_ref_channels.dropna(how='any')
        in_batch_overlaps[batch] = temp_ref_channels_filtered.shape[0] / temp_ref_channels.shape[0] * 100
        in_batch_counts[batch] = temp_ref_channels.count(axis=1)

    between_batches_overlap = ref_channels.dropna(how='any').shape[0]
    batch_counts = pd.DataFrame(in_batch_counts)
    # create protein annotation of occurrence percentage in ref channels in batches
    occ_in_percent_batches = (batch_counts[batch_counts == 0].isna().sum(axis=1) / batch_counts.shape[1]) * 100
    occ_in_percent_ref_channels = (batch_counts.sum(axis=1) / (batch_counts.shape[1] * 3)) * 100
    batch_counts = batch_counts.add_prefix('Batch')
    occ_in_percent_batches.name = 'Occurrence in % of batches'
    occ_in_percent_ref_channels.name = 'Occurrence in % of channels'
    intensity_range[0].name, intensity_range[
        1].name, median_intensity.name, std_intensity.name = 'Minimum intensity', 'Maximum intensity', 'Median intensity', 'Std intensity'
    combined_df = pd.concat(
        [batch_counts, occ_in_percent_batches, occ_in_percent_ref_channels, intensity_range[0], intensity_range[1],
         median_intensity, std_intensity, annot_channels],
        axis=1)
    combined_df.to_csv(os.path.join(results_folder, f'extracted_ref_channel_info_{data_type}.csv'))
    return combined_df


def get_ref_channels(df: pd.DataFrame):
    df = df.loc[:, df.columns.str.startswith('Reporter intensity')]
    return df


def get_annot(df: pd.DataFrame):
    df = df.loc[:, ['basket', 'basket_weights', 'sub_basket', 'sub_basket_weights']]
    return df


def get_intensity_channels(df: pd.DataFrame):
    df = df.loc[:, ~df.columns.str.startswith('Reporter intensity')]
    df = df.drop(columns=['basket', 'basket_weights', 'sub_basket', 'sub_basket_weights'])
    return df


def create_qc_report(data_dict: Dict[str, pd.DataFrame],
                     meta_data_dicts: Dict[str, pd.DataFrame]) -> None:
    # Get number of identifications + average Z-scores
    median_vs_count_dfs = get_z_score_idents([data_dict['full_proteome_measures_z'], data_dict['phospho_measures_z']],
                                             meta_data_dicts[
                                                 'Batch_No'])
    # pct_threshold = [1, 0.5, 0.1, 0.05] # slow
    pct_threshold = [1, 0.5]
    pca_colors, pca_markers = get_pca_color_marker(data_dict['preprocessed_fp'], meta_data_dicts['Batch_No'])
    pca_objects_list, _ = get_pca_objects([data_dict['preprocessed_fp'], data_dict['preprocessed_pp']],
                                          pct_threshold=pct_threshold,
                                          do_pca=True)

    targets = get_pca_targets(data_dict['preprocessed_fp'], meta_data_dicts['Batch_No'])
    principal_dfs = get_pca_dfs([data_dict['preprocessed_fp'], data_dict['preprocessed_pp']], pca_objects_list, targets,
                                pca_type='both')
    pc_variances = [[], []]
    for i, data_type in enumerate(pca_objects_list):
        for j, pca_object in enumerate(data_type):
            if j == 0:
                pc_variances[i].append(pca_object.explained_variance_ratio_)
            else:
                pc_variances[i].append(pca_object.var_exp)

    # TODO in plot function and principal df general change 'Batch' to label or so
    # TODO in plot function - title dependent on which label
    plot_qc_results(results_folder, median_vs_count_dfs, principal_dfs, pc_variances, pca_colors, pca_markers)


def investigate_pca_feature_importance(fp_df: pd.DataFrame,
                                       pca_objects_list: List[Any]) -> pd.DataFrame:
    """Retrieve per protein impact on PC and save to table. Two methods used"""
    # PPCA object information
    fp_df = fp_df.loc[fp_df.count(axis=1) >= fp_df.shape[1] * 0.5, :]
    # feature_importance = pd.DataFrame(pca_objects_list[0].C, columns=['PC1', 'PC2'])
    feature_importance = pd.DataFrame(pca_objects_list.C, columns=['PC1', 'PC2'])
    feature_importance['PC1_abs'] = abs(feature_importance['PC1'])
    feature_importance['PC2_abs'] = abs(feature_importance['PC2'])
    feature_importance.index = fp_df.index
    feature_importance_pc = feature_importance.sort_values('PC1_abs', ascending=False)
    return feature_importance_pc


def read_raw(file: Union[str, Path],
             sample_annotation_df: pd.DataFrame,
             data_type: str = 'fp') -> pd.DataFrame:
    """Read raw preprocessing files for investigation (before picked group FDR)"""
    df = pd.read_csv(file)
    if data_type == 'fp':
        df = df.set_index(['Batch', 'Gene names', 'Proteins'])
    else:
        df = df.set_index(['Batch', 'Gene names', 'Modified sequence'])
    df = df[df['Potential contaminant'] != '+']
    df = df[df['Reverse'] != '+']
    df = df.filter(regex='Reporter intensity corrected \d{1,2}')
    df = df.reset_index()
    df['new_genes'] = df['Gene names']
    df = prep.sum_peptide_intensities(df, debug=True)
    df = prep.convert_long_to_wide_format(df, debug=True)
    channel_to_sample_id_dict = bin.sample_annotation.get_channel_to_sample_id_dict(sample_annotation_df, remove_qc_failed=False)
    df = prep.rename_columns_with_sample_ids(df, channel_to_sample_id_dict, index_cols=['Gene names', 'Proteins'])
    columns = [column.strip() for column in df.columns]
    df.columns = columns
    return df


def get_batch_from_channel_name(x: str,
                                must_contain_batch: bool,
                                ref: bool = False) -> str:
    """Retrieve Batch identifier in channel name"""
    # elif must_contain_batch:
    #     match = re.search(r'\d{1,2} Batch\d{1,2}', x)
    #     x = match.group(0)
    if not must_contain_batch:
        if 'Batch' in x:
            match = re.search(r'Batch\d{1,2}', x)
            x = match.group(0)
        else:
            x = ' '
    if must_contain_batch:
        if 'Batch' in x:
            match = re.search(r'\d{1,2} [a-zA-Z]+_Batch\d{1,2}', x)
            x = match.group(0)
        else:
            if ref:
                x = ' '
            else:
                x = x
    return x


def get_unique_ordered_list(seq) -> List:
    """Return ordered set"""
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]


def get_z_score_idents(dfs: List[pd.DataFrame],
                       meta_data_df: pd.DataFrame) -> List[pd.DataFrame]:
    """Return DataFrame with median z-score and # identifications per patient"""
    median_vs_count_dfs = []
    for i, df in enumerate(dfs):
        # describe
        df = df.filter(regex=r'zscore_', axis=1)
        sample_names = df.columns.str.lstrip('zscore_')
        df = df.loc[:, sample_names.isin(meta_data_df['Sample name'].values.tolist())]

        median_vs_count_df = pd.DataFrame(pd.concat([df.median(axis=0), df.count()], axis=1))
        # remove rows not samples
        median_vs_count_df = median_vs_count_df.rename(columns={0: 'Median z-score', 1: 'Number of IDs'})
        samples = pd.Series(median_vs_count_df.index)
        samples = [sample.split('_')[1].strip() for sample in samples]
        median_vs_count_df.index = samples
        batch = [meta_data_df.loc[meta_data_df['Sample name'] == sample, 'Batch Name'].values.tolist()[0] for sample in
                 median_vs_count_df.index]
        median_vs_count_df['Batch'] = batch
        median_vs_count_dfs.append(median_vs_count_df)
    return median_vs_count_dfs


def plot_qc_results(results_folder: Union[str, Path],
                    median_vs_count_dfs: List[pd.DataFrame],
                    principal_dict: List[pd.DataFrame],
                    principal_var_dict: List,
                    new_batch: int) -> None:
    """QC report with plot of Median Z-score vs number of identifications, Proteins/p-sites per sample, Cohort PCA plot w different
    thresholds ...  """

    with PdfPages(os.path.join(results_folder, 'Results_investigation_plots', 'QC_report.pdf')) as pdf:

        # Generate title page
        firstPage = plt.figure(figsize=(8.27, 11.69))
        firstPage.clf()
        txt = 'QC plots of sample processing results'
        txt2 = 'for TOPAS WP3 pipeline'
        firstPage.text(0.5, 0.5, txt, transform=firstPage.transFigure, size=24, ha="center")
        firstPage.text(0.5, 0.45, txt2, transform=firstPage.transFigure, size=12, ha="center")
        pdf.savefig()
        plt.close()

        # Prepare colors (+marker symbols) for different plots
        fp, pp = principal_dict['fp'], principal_dict['pp']
        fp_pc_variances, pp_pc_variances = principal_var_dict['fp'], principal_var_dict['pp']
        # get color symbols
        colors_symbols = pca_tools.get_pca_colors_symbols(fp[0]['Batch_No'])
        targets = median_vs_count_dfs[0]['Batch'].unique()
        # TODO: move these plots to functions

        # Median Z-score vs number of identifications
        for i, df in enumerate(median_vs_count_dfs):
            if i == 0:
                fig, ax = plt.subplots(2, 1, figsize=(8.27, 11.69))
                for target, marker_tuple in zip(targets, colors_symbols):
                    indicesToKeep = df['Batch'] == target

                    ax[0].scatter(x=df.loc[indicesToKeep, 'Median z-score'], y=df.loc[indicesToKeep, 'Number of IDs'],
                                  c=marker_tuple[0],
                                  marker=marker_tuple[1],
                                  s=30,
                                  alpha=0.9)
                r1, _ = stats.pearsonr(df['Median z-score'], df['Number of IDs'])
            else:
                for target, marker_tuple in zip(targets, colors_symbols):
                    indicesToKeep = df['Batch'] == target
                    ax[1].scatter(x=df.loc[indicesToKeep, 'Median z-score'], y=df.loc[indicesToKeep, 'Number of IDs'],
                                  c=marker_tuple[0],
                                  marker=marker_tuple[1],
                                  s=30,
                                  alpha=0.9)
                r2, _ = stats.pearsonr(df['Median z-score'], df['Number of IDs'])
                ax[0].vlines(0, ax[0].get_ylim()[0], ax[0].get_ylim()[1], colors='k', linestyles='dashed')
                ax[0].vlines(-0.2, ax[0].get_ylim()[0], ax[0].get_ylim()[1], colors='grey', linestyles='dashed')
                ax[0].vlines(0.2, ax[0].get_ylim()[0], ax[0].get_ylim()[1], colors='grey', linestyles='dashed')
                ax[1].vlines(0, ax[1].get_ylim()[0], ax[1].get_ylim()[1], colors='k', linestyles='dashed')
                ax[1].vlines(-0.2, ax[1].get_ylim()[0], ax[1].get_ylim()[1], colors='grey', linestyles='dashed')
                ax[1].vlines(0.2, ax[1].get_ylim()[0], ax[1].get_ylim()[1], colors='grey', linestyles='dashed')
                x_corner1 = ax[0].get_xlim()[1] * 0.5
                y_corner1 = ax[0].get_ylim()[1] * 0.05
                x_corner2 = ax[1].get_xlim()[1] * 0.05
                y_corner2 = ax[1].get_ylim()[1] * 0.5
                ax[0].text(ax[0].get_xlim()[1] - x_corner1, ax[0].get_ylim()[1] - y_corner1, '$\\rho$=' + str(round(r1, 2)))
                ax[1].text(ax[1].get_xlim()[1] - x_corner2, ax[1].get_ylim()[1] - y_corner2, '$\\rho$=' + str(round(r2, 2)))
                ax[0].legend(df['Batch'].unique(), prop={'size': 5})
                ax[1].legend(df['Batch'].unique(), prop={'size': 5})
                plt.xlim([-1.25, 1.25])
                plt.title('Correlation between Z-score and # identifications (top FP, bottom PP)')
                plt.tight_layout()
                pdf.savefig(fig)
                plt.close()

        # Proteins/p-sites per sample
        for i, df in enumerate(median_vs_count_dfs):
            # TODO small functions for each plot... same structure with each batch one color/marker
            # TODO make input decide if highlight one batch
            df = df.reset_index()
            colors_batch_dict = {}
            # TODO: make sure this works no matter if string or int.. if string expected and int given make it a string
            for batch in df['Batch']:
                if int(batch) == new_batch:
                    colors_batch_dict[batch] = 'red'
                else:
                    colors_batch_dict[batch] = 'silver'
            if i == 0:
                fig, axes = plt.subplots(2, 1, figsize=(8.27, 11.69))
                axes[0].bar(x=df['index'], height=df['Number of IDs'], color=df['Batch'].map(colors_batch_dict), alpha=0.9)
            else:
                axes[1].bar(x=df['index'], height=df['Number of IDs'], color=df['Batch'].map(colors_batch_dict), alpha=0.9)
                plt.title('# Proteins/peptides per sample (top FP proteins, bottom PP peptides)')
                axes[0].tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
                axes[1].tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
                plt.tight_layout()
                pdf.savefig(fig)
                plt.close()

        # Reference channels ID numbers
        # should be a barplot? reuse code from above

        # Reference channels quant stability


        # PCA plots
        description = ['PCA plot of cohort subsetted to no missing values', 'PPCA plot of cohort subsetted to no missing value',
                       'PPCA with threshold >50% of patients', 'PPCA with threshold >10% of patients',
                       'PPCA with threshold >5% of patients']
        title = 'PCA Batches (top FP, bottom PP)'
        description = 'PPCA with threshold >50% of patients'
        # Loop through PCA plots making two subplots (FP, PP)
        for i in range(len(fp)):
            fp_pc_variances_temp = fp_pc_variances[i]
            pp_pc_variances_temp = pp_pc_variances[i]
            fp_principal_dfs_temp = fp[i]
            pp_principal_dfs_temp = pp[i]
            fig, ax = plt.subplots(2, 1, figsize=(8.27, 11.69))
            fig = prepare_double_plot(fig, ax, df=[fp_principal_dfs_temp, pp_principal_dfs_temp],
                                      exp_var=[fp_pc_variances_temp, pp_pc_variances_temp], target='Batch_No',
                                      color_symbols=colors_symbols,
                                      suptitle=description, title=title)
            pdf.savefig(fig)
            plt.close()
            # TODO check if it works with multiple input dfs like different pct

        # Create PCA plot with new batch samples highlighted on grey background cohort
        colors_symbols = [(color, 'o') for color in pd.Series(df['Batch'].unique()).map(colors_batch_dict)]
        fig, ax = plt.subplots(2, 1, figsize=(8.27, 11.69))
        fig = prepare_double_plot(fig, ax, df=[fp_principal_dfs_temp, pp_principal_dfs_temp],
                                  exp_var=[fp_pc_variances_temp, pp_pc_variances_temp], target='Batch_No',
                                  color_symbols=colors_symbols,
                                  suptitle=description, title=title)
        pdf.savefig(fig)
        plt.close()


def prepare_double_plot(fig,
                        ax,
                        df: List,
                        exp_var: List,
                        target: str,
                        color_symbols: List[Tuple],
                        suptitle: str = None,
                        title: str = None):
    targets = df[0][target].unique()

    for temp_target, marker_tuple in zip(targets, color_symbols):
        temp_df_fp = df[0][df[0][f'{target}'] == temp_target]
        temp_df_pp = df[1][df[1][f'{target}'] == temp_target]

        ax[0].scatter(temp_df_fp[['Principal component 1']],
                      temp_df_fp[['Principal component 2']], s=30, alpha=1, c=marker_tuple[0], marker=marker_tuple[1])
        ax[1].scatter(temp_df_pp[['Principal component 1']],
                      temp_df_pp[['Principal component 2']], s=30, alpha=1, c=marker_tuple[0], marker=marker_tuple[1])
        ax[0].set_xlabel('PC1 (%.0f%%)' % (np.round(exp_var[0][0] * 100)), fontsize=14)
        ax[0].set_ylabel('PC2 (%.0f%%)' % (np.round(exp_var[0][1] * 100)), fontsize=14)
        ax[1].set_xlabel('PC1 (%.0f%%)' % (np.round(exp_var[1][0] * 100)), fontsize=14)
        ax[1].set_ylabel('PC2 (%.0f%%)' % (np.round(exp_var[1][1] * 100)), fontsize=14)
        # ax[0].tick_params(axis='x', labelsize=13)
        # ax[1].tick_params(axis='x', labelsize=13)
        # ax[0].tick_params(axis='y', labelsize=13)
        # ax[1].tick_params(axis='y', labelsize=13)

    if title:
        plt.title(title, fontsize=12)
    if suptitle:
        plt.suptitle(suptitle, fontsize=16)
    plt.tight_layout()
    return fig


def read_patient_data(folder: Union[str, Path], sample_annotation_df: pd.DataFrame, data_types: List[str]) -> Dict[
    str, pd.DataFrame]:
    """
    Read data files of kind 'preprocessed_<data_type>_ref.csv', 'preprocessed_<data_type>.csv', 'debug_preprocess_<data_type>_raw.csv',
    'debug_preprocess_<data_type>_before_ms1.csv', 'debug_preprocess_<data_type>_after_ms1.csv', 'basket_scores_long.tsv',
    'full_proteome_measures_z.tsv'.......   <data_type> being either fp or pp
    :param folder: folder to find patient data in
    :param sample_annotation_df: annotation to get sample names mapping for debug files
    :return: dictionary with prepared dataframes using filename without ending as keys
    """
    data_dict = defaultdict()
    # results_files = ['debug_preprocessed_pp2_complete_raw.csv', 'debug_preprocessed_pp2_before_ms1_correction.csv',
    #                  'debug_preprocessed_pp2_after_ms1_correction.csv']
    # results_files = ['preprocessed_fp_ref.csv', 'preprocessed_pp_ref.csv', 'preprocessed_fp.csv', 'preprocessed_pp.csv',
    #                  'basket_scores_long.tsv', 'full_proteome_measures_rank.tsv', 'phospho_measures_rank.tsv',
    #                  'full_proteome_measures_z.tsv', 'phospho_measures_z.tsv']
    # , 'debug_preprocessed_fp2_complete_raw.csv',
    #            'debug_preprocessed_fp2_before_ms1_correction.csv', 'debug_preprocessed_fp2_after_ms1_correction.csv']
    # results_files = [f'preprocessed_{data_type}_with_ref.csv' for data_type in data_types] + \
    #                 [f'preprocessed_{data_type}.csv' for data_type in data_types] + \
    #                 [f'annot_{data_type}.csv' for data_type in data_types] + \
    #                 [f'annot_{data_type}_with_ref.csv' for data_type in data_types]

    results_files = [f'annot_{data_type}_with_ref.csv' for data_type in data_types]

    # z_score_measures = ['full_proteome_measures_z.tsv', 'phospho_measures_z.tsv']
    # for measure in z_score_measures:
    #     results_files.append(measure)

    # Look through files and if they exist add them to dict
    for result_file in results_files:
        logger.info(f'Reading in: {result_file}')
        print(f'Reading in: {result_file}')
        if os.path.exists(os.path.join(folder, result_file)):
            modality = result_file.split('.')[0]

            # if debug in it run prep first
            if 'debug' in result_file:
                if 'fp' in result_file:
                    df = read_raw(os.path.join(folder, result_file), sample_annotation_df, data_type='fp')
                if 'pp' in result_file:  # requires high memory and possibly CPU
                    df = read_raw(os.path.join(folder, result_file), sample_annotation_df, data_type='pp')
            elif 'basket' in result_file:
                df = pd.read_csv(os.path.join(folder, result_file), sep='\t', index_col='Sample')
            elif '.tsv' in result_file:
                df = pd.read_csv(os.path.join(folder, result_file), sep='\t')
            else:
                df = pd.read_csv(os.path.join(folder, result_file))
            if 'fp' in result_file or 'full' in result_file:
                df = df.set_index('Gene names')
            elif 'pp' in result_file and 'preprocessed' or 'annot' in result_file:
                df = df.set_index(['Modified sequence', 'Gene names', 'Proteins'])
            elif 'phospho' in result_file:
                df = df.set_index('Modified sequence')
            if 'basket' not in result_file and 'ref' not in result_file:
                df = df.filter(regex=r'[A-Z,a-z]+\d{1,3}-\S+-\S')
            df = df[df.columns.drop(list(df.filter(regex='metadata')))]
            data_dict[modality] = df
    return data_dict


if __name__ == "__main__":
    # Parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", dest="results", default=None, metavar="DIR", required=True, help="Path to results to check")

    argv = sys.argv[1:]
    args = parser.parse_args(argv)
    results_folder = Path(args.results)

    config_file = os.path.join(results_folder, 'configs.json')
    configs = config.load(config_file)

    check_results(results_folder, configs['new_batch'], configs['sample_annotation'], configs['metadata_annotation'],
                  configs["data_types"])
