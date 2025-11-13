import os
import sys
from pathlib import Path
import argparse
import logging

import numpy as np
import pandas as pd

from .. import config
from .. import utils
from . import scoring
from topas_pipeline import phospho_grouping
from topas_pipeline import bridge_normalization
from . import rtk_substrate_phosphorylation as rtk_scoring
from . import ck_substrate_phosphorylation as ck_scoring

# hacky way to get the package logger instead of just __main__ when running as a module
logger = logging.getLogger(__package__ + "." + Path(__file__).stem)


def protein_phospho_scoring_peptidoforms(
    results_folder: str, sample_annotation_file: str, metadata_file: str
):
    logger.info("Running protein phosphorylation scoring module")
    results_folder = Path(results_folder)
    protein_phosphorylation_file = (
        results_folder / "topas_scores" / "protein_phosphorylation_scores.tsv"
    )
    if protein_phosphorylation_file.is_file():
        logger.info(
            f"Protein phosphorylation scoring skipped - found files already processed"
        )
        return

    cohort_intensities_df = phospho_grouping.read_cohort_intensities_df(
        results_folder, sample_annotation_file
    )

    cohort_batch_corrected_df = bridge_normalization.read_cohort_batch_corrected_df(
        results_folder
    )

    cohort_intensities_df = rtk_scoring.update_with_batch_corrected_intensities(
        cohort_intensities_df, cohort_batch_corrected_df
    )

    annotation_df = cohort_intensities_df.index.to_frame().rename(
        columns={"Gene names": "Phosphoprotein"}
    )["Phosphoprotein"]
    protein_phosphorylation_score_df = (
        ck_scoring.compute_substrate_phosphorylation_scores(
            cohort_intensities_df, annotation_df, explode=False
        )
    )

    ck_scoring.save_scores(
        protein_phosphorylation_score_df, protein_phosphorylation_file
    )
    if metadata_file is not None and len(metadata_file) > 0:
        ck_scoring.save_scores(
            protein_phosphorylation_score_df, protein_phosphorylation_file, metadata_file
        )


def protein_phospho_scoring(results_folder, preprocessed_protein_df):
    logger.info("Running protein phosphorylation scoring module")

    if os.path.exists(results_folder / "protein_results" / "protein_scores.tsv"):
        logger.info(
            f"Protein phosphorylation scoring skipped - found files already processed"
        )
        return

    if not os.path.exists(protein_folder := results_folder / "protein_results"):
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
        results_folder / "protein_results" / "protein_scores_transposed.tsv",
        sep="\t",
        float_format="%.4g",
    )


@utils.validate_file_access
def load_protein_phosphorylation(
    results_folder,
    protein_results_folder: str = "topas_scores",
    remove_multi_gene_groups: bool = False,
):
    protein_scores_file = (
        results_folder / protein_results_folder / "protein_phosphorylation_scores.tsv"
    )
    protein_phosphorylation_df = pd.read_csv(protein_scores_file, index_col=0, sep="\t")
    protein_phosphorylation_df = protein_phosphorylation_df.T
    protein_phosphorylation_df.index.name = "Gene names"

    if remove_multi_gene_groups:
        # remove phosphoprotein groups which are a result of shared peptides
        protein_phosphorylation_df = protein_phosphorylation_df.loc[
            ~protein_phosphorylation_df.index.str.contains(";")
        ]
    return protein_phosphorylation_df


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", help="Absolute path to configuration file.")
    args = parser.parse_args(sys.argv[1:])
    configs = config.load(args.config)

    # preprocessed_df = scoring.topas_score_preprocess(configs.results_folder)
    # protein_phospho_scoring(configs.results_folder, preprocessed_df)

    protein_phospho_scoring_peptidoforms(
        configs.results_folder, configs.sample_annotation, configs.metadata_annotation
    )
