import os
import sys
import re
import argparse
from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib
from scipy import stats

import bin.basket_scoring as basket_scoring

if sys.platform.startswith('linux'):
    matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams['text.usetex'] = False
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns

# pd.set_option('display.max_rows', None)


def compare_results(older_results, newer_results):

    # if same and new generation do comparison of most top scoring?
    # significant ones

    # create directory for comparison results
    comparison_str = f"comparison_{str(older_results).split('/')[-1]}_to_{str(newer_results).split('/')[-1]}"
    # TODO: make not hardcoded the output_folder
    output_folder = f"{newer_results}/basket_score_comparisons"
    os.makedirs(os.path.join(output_folder, comparison_str), mode=0o700, exist_ok=True)

    score_correlations_path = os.path.join(output_folder, comparison_str, 'score_correlations_plots.pdf')
    sub_basket_score_correlations_path = os.path.join(output_folder, comparison_str, 'subbasket_score_correlations_plots.pdf')
    measure_correlations_path = os.path.join(output_folder, comparison_str, 'measures_correlations_plots.pdf')
    # comparison_files = [score_correlations_path, measure_correlations_path]
    comparison_files = [sub_basket_score_correlations_path, measure_correlations_path]
    # comparison_files = [measure_correlations_path, score_correlations_path]

    for file in comparison_files: # we dont need a loop I think

        if 'measures' in file:  # TODO remove?
            continue
        #     # Get measures
        #     measures_dict_1 = read_measures(older_results)
        #     measures_dict_2 = read_measures(newer_results)
        #
        #     # Get differences
        #     measure_changes = calc_changes(measures_dict_1, measures_dict_2)
        #
        #     # make two distinct functions that gets it for all measures
        #     diff_values, outliers, measures = [], [], []
        #     for i, measure in enumerate(measure_changes.keys()):
        #         print('Processing: ', measure)
        #         measures.append(measure)
        #
        #         # Get differences for histogram
        #         diff_values.append(get_diffs(measure_changes[measure]))
        #         # TODO - add min and max to plot or description  ?
        #
        #         # Get outliers for table
        #         outliers.append(get_outliers(measure_changes[measure]))
        #
        #     plot_to_pdf(diff_values=diff_values, file=file, measures=measures, outliers_rank=outliers)
        #
        #         # if 'z_score' in measure:
        #         #
        #         #     # not using inf and nan values to calculate mean and std
        #         #     mean_diff = np.ma.masked_invalid(diff_values).mean()
        #         #     std_diff = np.ma.masked_invalid(diff_values).std()
        #         #     # print(min_diff, max_diff)
        #         #     print(pd.Series(np.ma.masked_invalid(diff_values)).sort_values(ascending=False).head(10))
        #         #     print(mean_diff)
        #             # print(std_diff)
        #         #     diff_values = pd.Series(diff_values)
        #         #     new_zscore = diff_values.copy()
        #         #     new_zscore = (new_zscore - mean_diff) / std_diff

        else:
            # TODO: First optimize outlier definition
            # TODO: Get a summary for page 1
            # TODO: Make hyperlinks from page 1 to according page
            if 'subbasket' in file:
                basket_scores_1 = basket_scoring.read_sub_basket_scores(older_results).transpose()
                basket_scores_1.index.name = 'Sample'
                basket_scores_2 = basket_scoring.read_sub_basket_scores(newer_results).transpose()
                basket_scores_2.index.name = 'Sample'

            else:
                basket_scores_1 = basket_scoring.read_basket_scores(older_results)
                basket_scores_2 = basket_scoring.read_basket_scores(newer_results)

            print(basket_scores_1.head())
            print(basket_scores_2.head())
            merged = pd.merge(left=basket_scores_1, right=basket_scores_2, left_on='Sample', right_on='Sample', how='inner')
            common_columns = [c for c in merged.columns if "_x" in c and c.replace("_x", "_y") in merged.columns]

            measure_changes = calc_changes_baskets(merged, common_columns)
            diff_values, outliers = [], []
            for i, measure in enumerate(measure_changes.keys()):

                # get list
                temp_outlier = get_outliers(measure_changes[measure], basket=True)
                # print(temp_outlier)
                outliers.append(merged.loc[temp_outlier, [measure, measure.replace("_x", "_y")]])
            # print(outliers[0])

            plot_to_pdf(diff_values=merged, file=file, measures=common_columns, outliers_baskets=outliers)


def get_diffs(diff_df, basket=False):
    diff_values = diff_df.values
    diff_values = [value for sublist in diff_values for value in sublist]
    diff_values = pd.Series(np.ma.masked_invalid(diff_values))
    return diff_values


def get_outliers(diff_df, basket=False):

    if basket:
        outlier_table = diff_df.sort_values().tail(10).index.values.tolist()
    else:
        diff_df = diff_df.reset_index()
        if 'Gene names' in diff_df.columns:
            diff_df = diff_df.melt(id_vars='Gene names', value_vars=diff_df.columns.tolist())
        else:
            diff_df = diff_df.melt(id_vars='Modified sequence', value_vars=diff_df.columns.tolist())

        outlier_table = pd.concat([diff_df.nlargest(5, 'value'), diff_df.nsmallest(5, 'value')], axis=0)
    return outlier_table


def plot_to_pdf(diff_values, file, measures=None, outliers_rank=None, outliers_baskets=None):
    with PdfPages(f'{file}') as pdf:
        firstPage = plt.figure(figsize=(8.27, 11.69))
        firstPage.clf()
        txt = 'Comparison between runs'
        txt2 = 'for TOPAS WP3 sample processing pipeline'
        firstPage.text(0.5, 0.5, txt, transform=firstPage.transFigure, size=24, ha="center")
        firstPage.text(0.5, 0.45, txt2, transform=firstPage.transFigure, size=12, ha="center")
        pdf.savefig()
        plt.close()

        if 'score_correlations' in file:

            for i, c in enumerate(measures):
                outlier_temp = outliers_baskets[i]
                outlier_temp = outlier_temp.reset_index()
                outlier_temp = outlier_temp.round(2)

                # Calc statistics on values not being nan
                nas = np.logical_or(diff_values[c].isnull(), diff_values[c.replace("_x", "_y")].isnull())
                r, p_value = stats.pearsonr(diff_values[c][~nas], diff_values[c.replace("_x", "_y")][~nas])

                if i % 2 == 0:
                    fig, ax = plt.subplots(2, 2, figsize=(8.27, 11.69))
                    sns.scatterplot(data=diff_values, x=c, y=c.replace("_x", "_y"), ax=ax[i % 2, i % 2])  # [0,0]
                    sns.scatterplot(data=outlier_temp, x=c, y=c.replace("_x", "_y"), ax=ax[i % 2, i % 2], color='r')
                    xlim_min, xlim_max = ax[0, 0].get_xlim()[0], ax[0, 0].get_xlim()[1]
                    x = np.arange(xlim_min, xlim_max, step=(xlim_max-xlim_min)/100)
                    ax[0, 0].plot(x, x, color='grey', linestyle='dashed')

                    # Outlier table
                    ax[i % 2, i % 2+1].table(cellText=np.vstack([outlier_temp.columns, outlier_temp.values]), bbox=[0, 0, 1, 1])
                    ax[i % 2, i % 2 + 1].set_xticks([])
                    ax[i % 2, i % 2 + 1].set_yticks([])

                    # Pearson correlation
                    ax[0, 0].text(0.05, 0.95, '$\\rho$=' + str(round(r, 2)), transform=ax[0, 0].transAxes, va='top')

                if i % 2 == 1:
                    sns.scatterplot(data=diff_values, x=c, y=c.replace("_x", "_y"), ax=ax[i % 2, i % 2-1])  # [1,0]
                    sns.scatterplot(data=outlier_temp, x=c, y=c.replace("_x", "_y"), ax=ax[i % 2, i % 2-1], color='r')
                    xlim_min, xlim_max = ax[1, 0].get_xlim()[0], ax[1, 0].get_xlim()[1]
                    x = np.arange(xlim_min, xlim_max, step=(xlim_max-xlim_min)/100)
                    ax[1, 0].plot(x, x, color='grey', linestyle='dashed')

                    # Outlier table
                    ax[i % 2, i % 2].table(cellText=np.vstack([outlier_temp.columns, outlier_temp.values]), bbox=[0, 0, 1, 1])
                    ax[i % 2, i % 2].set_xticks([])
                    ax[i % 2, i % 2].set_yticks([])

                    # Pearson correlation
                    ax[1, 0].text(0.05, 0.95, '$\\rho$=' + str(round(r, 2)), transform=ax[1, 0].transAxes, va='top')

                    plt.tight_layout()
                    pdf.savefig(fig)
                    plt.close()

        if 'measures_correlations' in file:
            # Rank differences
            fig, ax = plt.subplots(2, 1, figsize=(8.27, 11.69))
            for i, vals in enumerate(diff_values[:2]):
                p = sns.histplot(data=vals, bins=50, ax=ax[i])
                p.set(yscale='log')
                ax[i].set_xlabel(f'Difference in {measures[:2][i]}')
                plt.tight_layout()
            pdf.savefig(fig)
            plt.close()
            # Rank outliers
            fig, ax = plt.subplots(2, 1, figsize=(8.27, 11.69))
            for i, vals in enumerate(outliers_rank[:2]):
                ax[i].axis('off')
                c = outliers_rank[:2][i].shape[1]
                ax[i].table(cellText=np.vstack([outliers_rank[i].columns, outliers_rank[i].values]), bbox=[0,0,1,1])
            pdf.savefig(fig)
            plt.close()
            # Z-score differences
            fig, ax = plt.subplots(2, 1, figsize=(8.27, 11.69))
            for i, vals in enumerate(diff_values[2:]):
                p = sns.histplot(data=vals, bins=50, ax=ax[i])
                p.set(yscale='log')
                ax[i].set_xlabel(f'Difference in {measures[2:][i]}')
                plt.tight_layout()
            pdf.savefig(fig)
            plt.close()
            # Z-score outliers
            fig, ax = plt.subplots(2, 1, figsize=(8.27, 11.69))
            for i, vals in enumerate(outliers_rank[2:]):
                ax[i].axis('off')
                c = outliers_rank[2:][i].shape[1]
                ax[i].table(cellText=np.vstack([outliers_rank[i].columns, outliers_rank[2:][i].values]), bbox=[0,0,1,1])
            pdf.savefig(fig)
            plt.close()


def calc_changes(dict1, dict2):
    difference_dict = {}
    for measure in dict1.keys():
        # find intersect
        intersect = np.intersect1d(dict1[measure].index, dict2[measure].index)
        dict1[measure] = dict1[measure].loc[intersect, :]
        dict2[measure] = dict2[measure].loc[intersect, :]
        difference_dict[measure] = dict1[measure]-dict2[measure]   # TODO - do we want abs(diff) ??
    return difference_dict


def calc_changes_baskets(df, common_columns):
    difference_dict = {}
    for c in common_columns:
        difference_dict[c] = abs(df[c] - df[c.replace("_x", "_y")])
    return difference_dict


def read_measures(folder):
    ranks_fp = pd.read_csv(os.path.join(folder, 'full_proteome_measures_rank.tsv'), sep='\t', index_col='Gene names')
    ranks_pp = pd.read_csv(os.path.join(folder, 'phospho_measures_rank.tsv'), sep='\t', index_col='Modified sequence')
    z_score_fp = pd.read_csv(os.path.join(folder, 'full_proteome_measures_z.tsv'), sep='\t', index_col='Gene names')
    z_score_pp = pd.read_csv(os.path.join(folder, 'phospho_measures_z.tsv'), sep='\t', index_col='Modified sequence')
    ranks_fp = ranks_fp.filter(regex=r'rank_[A-Z]\d+')
    ranks_pp = ranks_pp.filter(regex=r'rank_[A-Z]\d+')
    z_score_fp = z_score_fp.filter(regex=r'zscore')
    z_score_pp = z_score_pp.filter(regex=r'zscore')
    return {'rank_fp': ranks_fp, 'rank_pp': ranks_pp, 'z_score_fp': z_score_fp, 'z_score_pp': z_score_pp}


# def read_basket_scores(folder):
#     basket_scores = pd.read_csv(os.path.join(folder, 'basket_scores.tsv'), sep='\t', index_col='Sample')
#     basket_scores = basket_scores.rename(columns=lambda x: x.lower())
#     return basket_scores


if __name__ == "__main__":
    # Parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--folder_1", dest="folder_1", default=None, metavar="DIR", required=True, help="Path to old results to use for "
                                                                                                        "comparison")
    parser.add_argument("--folder_2", dest="folder_2", default=None, metavar="DIR", required=True, help="Path to new results to use for "
                                                                                                        "comparison")

    argv = sys.argv[1:]
    args = parser.parse_args(argv)
    folder_1, folder_2 = Path(args.folder_1), Path(args.folder_2)

    # TODO error testing this also that it takes 2!! arguments
    # Ensure a config was passed to the script
    # if not args.folder_1:
    #     print("No folders provided.")
    #     exit()
    # else:

    compare_results(folder_1, folder_2)
