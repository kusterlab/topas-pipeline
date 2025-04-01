import os
import sys
from pathlib import Path
import argparse
import logging

import numpy as np
import pandas as pd

from .. import config
from . import scoring

# hacky way to get the package logger instead of just __main__ when running as a module
logger = logging.getLogger(__package__ + "." + Path(__file__).stem)


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

    protein_df = scoring.cap_zscores_and_weights(preprocessed_protein_df)

    logger.info("  Calculate protein phosphorylation scores")
    score_dataframe = protein_df.groupby("Gene names")
    score_dataframe = score_dataframe.agg(
        **(
            {
                score: pd.NamedAgg(column=score, aggfunc="sum")
                for score in protein_df.columns
                if "pat_" in score
            }
        )
    ).reset_index()
    score_dataframe = score_dataframe.replace(0, np.nan)

    logger.info("  2nd level z-scoring, adding target space and writing results")

    # what is happening here?
    protein_scores = scoring.second_level_z_scoring(score_dataframe, "Gene names")
    protein_scores = score_dataframe
    protein_scores.to_csv(
        os.path.join(results_folder, "protein_results", "protein_scores.tsv"),
        sep="\t",
        index=False,
        float_format="%.4g",
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
        float_format="%.4g",
    )


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", help="Absolute path to configuration file.")
    args = parser.parse_args(sys.argv[1:])
    configs = config.load(args.config)

    preprocessed_df = scoring.topas_score_preprocess(configs.results_folder)
    protein_phospho_scoring(configs.results_folder, preprocessed_df)
