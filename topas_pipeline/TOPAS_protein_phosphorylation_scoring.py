# %%
import os
import sys
from pathlib import Path
import argparse
import logging

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import topas_pipeline.config as config
import topas_pipeline.TOPAS_scoring_functions as scoring

# hacky way to get the package logger instead of just __main__ when running as a module
logger = logging.getLogger(__package__ + "." + Path(__file__).stem)

pd.set_option("display.max_columns", None)
pd.set_option("display.width", None)


def read_protein_scoring(results_folder: str):
    if os.path.exists(
        os.path.join(results_folder, "protein_results", "protein_scores.tsv")
    ):
        try:
            protein_scores = os.path.join(
                results_folder, "protein_results", "protein_scores.tsv"
            )
            protein_scores_df = pd.read_csv(
                protein_scores, index_col="Gene names", sep="\t"
            )
            protein_scores_df = protein_scores_df.filter(regex="pat_")
        except PermissionError:
            raise PermissionError(
                f"Cannot open protein phosphorylation scores file, check if you have it open in Excel. {protein_scores}"
            )
    return protein_scores_df


def protein_score_preprocess(results_folder):
    filepath = os.path.join(results_folder, "protein_score_preprocessed.tsv")
    if os.path.exists(filepath):
        return pd.read_csv(filepath, sep="\t")

    patients_zscores = pd.read_csv(
        os.path.join(results_folder, "phospho_measures_z.tsv"),
        keep_default_na=False,
        sep="\t",
    )
    patients_zscores = patients_zscores.rename(
        columns={
            colname: colname.replace("zscore_", "pat_")
            for colname in patients_zscores.columns.tolist()
        }
    )
    patients_list = patients_zscores.filter(regex="pat_").columns.tolist()
    patients_zscores = patients_zscores.loc[
        :, ["Gene names", "Modified sequence", "Proteins"] + patients_list
    ]
    patients = patients_zscores

    patients.set_index("Modified sequence", inplace=True)
    patients = patients.dropna(subset="Gene names")
    patients = scoring.calculate_peptide_occurrence(patients)
    patients["Gene names"] = patients["Gene names"].str.split(";")
    patients = patients.explode("Gene names")
    patients.to_csv(filepath, sep="\t")
    return patients


def protein_phospho_scoring(results_folder, preprocessed_protein_df):
    logger.info("Running protein phosphorylation scoring module")

    if os.path.exists(
        os.path.join(results_folder, "protein_results", "protein_scores.tsv")
    ):
        logger.info(
            f"Protein phosphorylation scoring skipped - found files already processed"
        )
        return

    if not os.path.exists(
        protein_folder := os.path.join(results_folder, "protein_results")
    ):
        os.makedirs(protein_folder)

    # TODO: Why are you not using the functions of psite_scoring_functions you moron?
    protein_df = preprocessed_protein_df.copy()
    patcols = [col for col in protein_df.columns if "pat_" in col]
    protein_df[patcols] = protein_df[patcols].astype(float).clip(upper=4, lower=-4)

    logger.info("  Calculate protein phosphorylation scores")
    score_dataframe = protein_df.groupby("Gene names")
    score_dataframe = score_dataframe.agg(
        **(
            {
                score: pd.NamedAgg(column=score, aggfunc=sum)
                for score in protein_df.columns
                if "pat_" in score
            }
        )
    ).reset_index()
    score_dataframe = score_dataframe.replace(0, np.nan)

    # score_dataframe = score_dataframe.rename(
    #     columns={col: col.replace('pat_', '') for col in score_dataframe.columns})

    # plot_histograms_to_check_normality(score_dataframe)

    logger.info("  2nd level z-scoring, adding target space and writing results")

    # what is happening here?
    protein_scores = scoring.second_level_z_scoring(score_dataframe, "Gene names")
    protein_scores = score_dataframe
    protein_scores.to_csv(
        os.path.join(results_folder, "protein_results", "protein_scores.tsv"),
        sep="\t",
        index=False,
    )
    protein_scores_t = (
        protein_scores.sort_values(by="Gene names")
        .set_index("Gene names")
        .transpose()
        .sort_index()
    )
    protein_scores_t.to_csv(
        os.path.join(
            results_folder, "protein_results", "protein_scores_transposed.tsv"
        ),
        sep="\t",
    )


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", help="Absolute path to configuration file.")
    args = parser.parse_args(sys.argv[1:])
    configs = config.load(args.config)

    preprocessed_protein_df = protein_score_preprocess(configs["results_folder"])
    protein_phospho_scoring(configs["results_folder"], preprocessed_protein_df)
