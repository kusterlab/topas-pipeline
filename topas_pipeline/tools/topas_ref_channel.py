import os
import re
from typing import Dict
import subprocess
import shutil

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

sns.set(font="Arial")
font_size = 12
sns.set_context(
    "paper",
    rc={
        "font.size": font_size,
        "axes.titlesize": font_size,
        "axes.labelsize": font_size,
        "legend.fontsize": font_size,
        "xtick.labelsize": font_size,
        "ytick.labelsize": font_size,
    },
)
plt.rcParams["svg.fonttype"] = "none"
sns.set_style("ticks")

# results_folder = '/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_study/2024.08.14_CJ_new_paper_cohort'
results_folder = '/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_study/2025.06.04_CJ_pancancer315_replicates'
# output_folder = '/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/MT/2025.05.30_TOPAS_ref_channel'
output_folder = '/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/CJ/2025.06.12_TOPAS_ref_channel'
os.makedirs(output_folder, exist_ok=True)


import csv

file_path = f"{results_folder}/annot_pp.csv"
columns_of_interest = ['PSP Kinases', 'Site positions', 'Site positions identified (MQ)']

with open(file_path, 'r', newline='', encoding='utf-8') as f:
    reader = csv.reader(f)
    header = next(reader)  # read only the first line

# Get indices (0-based)
indices = {col: header.index(col) for col in columns_of_interest}

print(indices) ## wrong by one so do from terminal 
# head -n 1 annot_pp.csv | tr ',' '\n' | nl


metadata_df = pd.read_excel(f"{results_folder}/METADATA_315_newQC_replicates.xlsx")
sample_batch_dict = metadata_df.set_index("Sample name")["Batch_No"].to_dict()

kinase_scores_with_ref_file = f"{output_folder}/kinase_results_original/kinase_scores.tsv"
# print('here')
# if not os.path.isfile(kinase_scores_with_ref_file):
#     from topas_pipeline.topas import phosphorylation

#     # Extract the annotation columns from annot_pp.csv necessary for psite_scoring
#     # !cut -f 1,2,3,5150,5159,5148,2586 -d',' {results_folder}/annot_pp.csv > {output_folder}/annot_pp.csv

#     input_file = os.path.join(results_folder, "annot_pp.csv")
#     output_file = os.path.join(output_folder, "annot_pp.csv")

#     cmd = f"cut -d',' -f1,2,3,5151,5160,5160,5149 \"{input_file}\" > \"{output_file}\""  # shit hardcoded.....
#     subprocess.run(cmd, shell=True, check=True)

#     # 17 minutes
#     phosphorylation.psite_scoring(output_folder, extra_kinase_annot="", data_types=["fp", "pp"])


#     # Define paths
#     kr_old = os.path.join(output_folder, "kinase_results")
#     pr_old = os.path.join(output_folder, "protein_results")
#     kr_new = os.path.join(output_folder, "kinase_results_original")
#     pr_new = os.path.join(output_folder, "protein_results_original")

#     # Move old directories if they exist
#     if os.path.exists(kr_old):
#         shutil.move(kr_old, kr_new)
#     if os.path.exists(pr_old):
#         shutil.move(pr_old, pr_new)

#     # Create new empty directories
#     os.makedirs(kr_old, exist_ok=True)
#     os.makedirs(pr_old, exist_ok=True)
#     # !mv {output_folder}/kinase_results {output_folder}/kinase_results_original
#     # !mv {output_folder}/protein_results {output_folder}/protein_results_original
#     # !mkdir {output_folder}/kinase_results
#     # !mkdir {output_folder}/protein_results




# kinase_df = pd.read_csv(f"{output_folder}/kinase_results_original/kinase_scores.tsv", sep='\t', index_col=[0, 1])
# kinase_df

# def compute_ref_zscores_topas_subscore(kinase_df: pd.DataFrame):
#     mean = kinase_df.filter(regex=r'^pat_(?!ref_)').apply(np.nanmean, axis=1, result_type="reduce")
#     stdev = kinase_df.filter(regex=r'^pat_(?!ref_)').apply(np.nanstd, axis=1, result_type="reduce")
#     return kinase_df.subtract(mean, axis=0).divide(stdev, axis=0)

# kinase_with_ref_df = compute_ref_zscores_topas_subscore(kinase_df)
# kinase_with_ref_df.to_csv(f"{output_folder}/kinase_results/kinase_scores.tsv", sep='\t')
# kinase_with_ref_df


# protein_phosphorylation_df = pd.read_csv(f"{output_folder}/protein_results_original/protein_scores.tsv", sep='\t', index_col=0)
# protein_phosphorylation_with_ref_df = compute_ref_zscores_topas_subscore(protein_phosphorylation_df)
# protein_phosphorylation_with_ref_df.to_csv(f"{output_folder}/protein_results/protein_scores.tsv", sep='\t')
# protein_phosphorylation_with_ref_df

# src_file = os.path.join(results_folder, "annot_fp.csv")
# dst_file = os.path.join(output_folder, "annot_fp.csv")

# if not os.path.isfile(f"{output_folder}/basket_scores_4th_gen_zscored.tsv"):
#     from topas_pipeline.topas import topas

#     # Copy the file
#     shutil.copyfile(src_file, dst_file)

#     # 3 minutes
#     topas.compute_topas_scores(
#         output_folder,
#         metadata_file="/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_MTBs_Evaluation/METADATA_315_newQC_replicates.xlsx",
#         topas_annotation_file="/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_MTBs_Evaluation/TOPASscores_POI_AS_250307.xlsx",
#     )


def compute_ref_zscores_topas_subscore(kinase_df: pd.DataFrame):
    mean = kinase_df.filter(regex=r'^pat_(?!ref_)').apply(np.nanmean, axis=1, result_type="reduce")
    stdev = kinase_df.filter(regex=r'^pat_(?!ref_)').apply(np.nanstd, axis=1, result_type="reduce")
    return kinase_df.subtract(mean, axis=0).divide(stdev, axis=0)


# now from here....
topas_df = pd.read_csv(f"{output_folder}/basket_scores_4th_gen_zscored.tsv", sep='\t', index_col=0)
topas_df = topas_df.replace(0, np.nan)
topas_with_ref_df = compute_ref_zscores_topas_subscore(topas_df.T).T

def get_qc_lot(batch_number: float):
    if batch_number < 73:
        return 1
    elif batch_number < 116:
        return 2
    else:
        return 3

def add_reference_information(topas_with_ref_df: pd.DataFrame):
    topas_with_ref_df['is_reference'] = topas_with_ref_df.index.str.startswith('pat_ref_')
    topas_with_ref_df['batch'] = np.nan
    topas_with_ref_df['qc_lot'] = np.nan
    topas_with_ref_df.loc[topas_with_ref_df['is_reference'], 'batch'] = topas_with_ref_df.loc[topas_with_ref_df['is_reference']].index.str.split('batch').str[1].astype(int)
    topas_with_ref_df.loc[topas_with_ref_df['is_reference'], 'qc_lot'] = topas_with_ref_df.loc[topas_with_ref_df['is_reference'], 'batch'].apply(get_qc_lot)
    return topas_with_ref_df

topas_with_ref_df = add_reference_information(topas_with_ref_df)
topas_with_ref_df

drop_cols = ["is_reference", "batch", "qc_lot"]
topas_stdevs = pd.concat(
    [
        topas_with_ref_df.filter(regex=r"^pat_(?!ref_)", axis=0)
        .drop(columns=drop_cols)
        .std(),
        topas_with_ref_df[topas_with_ref_df['is_reference']].drop(columns=drop_cols).std(),
        # topas_with_ref_df[topas_with_ref_df['is_reference'] & (topas_with_ref_df['qc_lot'] == 1)].drop(columns=drop_cols).std(),
        # topas_with_ref_df[topas_with_ref_df['is_reference'] & (topas_with_ref_df['qc_lot'] == 2)].drop(columns=drop_cols).std(),
        topas_with_ref_df[topas_with_ref_df['is_reference'] & (topas_with_ref_df['qc_lot'] == 3)].drop(columns=drop_cols).std(),
    ],
    axis=1,
)
# topas_stdevs.columns = ["patient channels", "reference channels all", "reference channels QC lot 1", "reference channels QC lot 2", "reference channels QC lot 3"]
topas_stdevs.columns = ["patient channels", "reference channels all", "reference channels QC lot 3"]
# sns.scatterplot(topas_stdevs)
# plt.xticks(rotation=90)
# plt.xlabel("TOPAS scores")
# plt.ylabel("standard deviation of z-scores")
# plt.ylim([0, 1.1])
# plt.savefig(f"{output_folder}/topas_stdevs.svg", bbox_inches="tight")
plt.figure(figsize=(18,40))
for i, topas_rtk in enumerate(topas_with_ref_df.drop(columns=drop_cols).columns):
    if topas_rtk == "is_reference":
        continue
    plt.subplot(7, 3, i+1)
    # plt.figure(figsize=(6,6))
    plt.title(topas_rtk)
    sns.swarmplot(topas_with_ref_df[[topas_rtk, 'is_reference']], y=topas_rtk, x='is_reference', hue='is_reference', size=3)
    # plt.savefig(f"{output_folder}/TOPAS_ref_channel_swarmplots_{topas_rtk}.png")
# plt.savefig(f"{output_folder}/TOPAS_ref_channel_swarmplots_merged.png")
plt.savefig(f"{output_folder}/TOPAS_ref_channel_swarmplots_merged.svg", format='svg')
