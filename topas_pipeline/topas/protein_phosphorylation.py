import sys
from pathlib import Path
import argparse
import logging

import pandas as pd

from .. import config
from .. import utils
from ..preprocess import phospho_grouping
from . import ck_substrate_phosphorylation as ck_scoring

# hacky way to get the package logger instead of just __main__ when running as a module
logger = logging.getLogger(__package__ + "." + Path(__file__).stem)


def protein_phospho_scoring(results_folder: str, metadata_file: str):
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
        f"{results_folder}/preprocessed_pp.csv",
        keep_identification_metadata_columns=False,
    )

    annotation_df = cohort_intensities_df.index.to_frame().rename(
        columns={"Gene names": "Phosphoprotein"}
    )["Phosphoprotein"]
    annotation_df = annotation_df.replace("", pd.NA).dropna()  # remove empty gene names
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
            protein_phosphorylation_score_df,
            protein_phosphorylation_file,
            metadata_file,
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

    protein_phospho_scoring(configs.results_folder, configs.metadata_annotation)
