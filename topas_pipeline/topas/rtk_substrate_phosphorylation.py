# compute substrate phosphorylation scores (a.k.a. kinase activity) for receptor
# tyrosine kinases using the manually curated confident substrates by Annika Schneider and Firas Hamood

import sys
import logging
from pathlib import Path
import os

import pandas as pd
import psite_annotation as pa

from .. import utils
from . import ck_substrate_phosphorylation
from topas_pipeline import phospho_grouping
from topas_pipeline import bridge_normalization

# hacky way to get the package logger instead of just __main__ when running as a module
logger = logging.getLogger(__package__ + "." + Path(__file__).stem)


def calculate_rtk_scores(
    results_folder: str,
    metadata_file: str,
    extra_kinase_annot: str,
    sample_annotation_file: str,
    fasta_file: str,
):
    results_folder = Path(results_folder)
    kinase_score_file = (
        results_folder / "topas_scores" / "rtk_substrate_phosphorylation_scores.tsv"
    )
    if kinase_score_file.is_file():
        logger.info(
            f"Receptor tyrosine substrate phosphorylation scoring skipped - found files already processed"
        )
        return
    cohort_annotated_sites_df = get_annotated_modified_sequence_groups(
        results_folder,
        extra_kinase_annot,
        fasta_file,
    )

    annotated_sites_mapping = get_annotated_sites_mapping(cohort_annotated_sites_df)

    annotated_cohort_intensities_df = read_annotated_cohort_intensities_df(
        results_folder,
        sample_annotation_file,
        cohort_annotated_sites_df,
    )

    annotated_modified_sequence_groups = (
        annotated_cohort_intensities_df.index.get_level_values(
            "Modified sequence group"
        )
    )
    annotated_cohort_batch_corrected_df = read_annotated_cohort_batch_corrected_df(
        results_folder, annotated_modified_sequence_groups
    )

    annotated_cohort_intensities_df = update_with_batch_corrected_intensities(
        annotated_cohort_intensities_df, annotated_cohort_batch_corrected_df
    )

    substrate_file = (
        results_folder / "topas_scores" / "rtk_substrate_peptide_intensities.tsv"
    )
    ck_substrate_phosphorylation.write_substrate_peptides(
        annotated_cohort_intensities_df, annotated_sites_mapping, substrate_file
    )

    scores = ck_substrate_phosphorylation.compute_substrate_phosphorylation_scores(
        annotated_cohort_intensities_df,
        annotated_sites_mapping,
    )

    ck_substrate_phosphorylation.save_scores(scores, kinase_score_file)
    if metadata_file is not None and len(metadata_file) > 0:
        ck_substrate_phosphorylation.save_scores(scores, kinase_score_file, metadata_file)


def get_annotated_modified_sequence_groups(
    results_folder: str,
    extra_kinase_annot: str,
    fasta_file: str,
):
    logger.info("Reading cohort modified sequence groups")
    # loading only these three columns takes 20 sec instead of 4 minutes
    pp_df = pd.read_csv(
        results_folder / "preprocessed_pp2_agg.csv",
        usecols=["Gene names", "Modified sequence group", "Proteins"],
    )

    # keep 0-based index as 'index' column so we only need to read the annotated rows later
    pp_df = pp_df.reset_index()
    pp_df_exploded = explode_modified_sequence_group(pp_df)

    logger.info("Annotating modified sequence groups with phosphorylating RTK")
    # annotate with extra_kinase_annot manually annotated RTK substrates
    cohort_annotated_sites_df = pa.addPeptideAndPsitePositions(
        pp_df_exploded[
            [
                "Modified sequence",
                "Modified sequence group",
                "Gene names",
                "Proteins",
                "index",
            ]
        ],
        fastaFile=fasta_file,
    )
    cohort_annotated_sites_df = pa.addPSPKinaseSubstrateAnnotations(
        cohort_annotated_sites_df, extra_kinase_annot
    )
    return cohort_annotated_sites_df


def explode_modified_sequence_group(df: pd.DataFrame):
    df["Modified sequence"] = df["Modified sequence group"].str.split(";")
    df = df.explode("Modified sequence")
    return df


def get_annotated_sites_mapping(df_patients_sites: pd.DataFrame):
    annotated_sites_mapping = df_patients_sites.loc[
        df_patients_sites["PSP Kinases"] != "",
        ["Modified sequence group", "PSP Kinases", "Gene names"],
    ]
    annotated_sites_mapping = annotated_sites_mapping.set_index(
        ["Modified sequence group", "Gene names"]
    )["PSP Kinases"]
    return annotated_sites_mapping


def read_annotated_cohort_intensities_df(
    results_folder: str,
    sample_annotation_file: str,
    cohort_annotated_sites_df: pd.DataFrame,
):
    logger.info("Read cohort intensities for annotated sites")
    keep_rows = cohort_annotated_sites_df.loc[
        cohort_annotated_sites_df["PSP Kinases"] != "", "index"
    ]
    skip_rows = cohort_annotated_sites_df.loc[
        ~cohort_annotated_sites_df["index"].isin(keep_rows), "index"
    ]
    skip_rows += 1  # add 1 to the row numbers to account for the header line

    return phospho_grouping.read_cohort_intensities_df(
        results_folder, sample_annotation_file, skip_rows
    )


def read_annotated_cohort_batch_corrected_df(
    results_folder: str, annotated_modified_sequence_groups: pd.Series
):
    logger.info("Read batch corrected cohort intensities for annotated sites")
    # only contains 6 out of ~50 annotated TOPAS-RTK sites
    df_patients = pd.read_csv(
        results_folder / "preprocessed_pp2_agg_batchcorrected.csv",
        usecols=["Gene names", "Modified sequence group"],
    )
    df_patients = df_patients.reset_index()
    skiprows = df_patients.loc[
        ~df_patients["Modified sequence group"].isin(
            annotated_modified_sequence_groups
        ),
        "index",
    ]
    skiprows += 1  # add 1 to the row numbers to account for the header line

    return bridge_normalization.read_cohort_batch_corrected_df(results_folder, skiprows)


def update_with_batch_corrected_intensities(
    cohort_intensities_df: pd.DataFrame, cohort_batch_corrected_df: pd.DataFrame
) -> pd.DataFrame:
    return pd.concat(
        [
            cohort_intensities_df.loc[
                ~cohort_intensities_df.index.isin(cohort_batch_corrected_df.index)
            ],
            cohort_batch_corrected_df.loc[
                cohort_batch_corrected_df.index.isin(cohort_intensities_df.index)
            ],
        ]
    )


@utils.validate_file_access
def load_substrate_phosphorylation(
    results_folder,
    kinase_results_folder: str = "topas_scores",
):
    kinase_score_file = (
        results_folder
        / kinase_results_folder
        / "rtk_substrate_phosphorylation_scores.tsv"
    )

    kinase_scores_df = pd.read_csv(
        kinase_score_file,
        index_col=0,
        sep="\t"
    )
    kinase_scores_df = kinase_scores_df.T
    kinase_scores_df.index.name = "PSP Kinases"
    kinase_scores_df.index = kinase_scores_df.index.str.replace(
        "RTK-TOPAS", "TOPAS", regex=False
    )

    return kinase_scores_df


"""
python3 -m topas_pipeline.topas.rtk_scoring -c config_patients.json
"""
if __name__ == "__main__":
    import argparse

    from .. import config

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-c", "--config", required=True, help="Absolute path to configuration file."
    )
    args = parser.parse_args(sys.argv[1:])

    configs = config.load(args.config)

    calculate_rtk_scores(
        results_folder=configs.results_folder,
        metadata_file=configs.metadata_annotation,
        extra_kinase_annot=configs.clinic_proc.extra_kinase_annot,
        sample_annotation_file=configs.sample_annotation,
        fasta_file=configs.preprocessing.fasta_file,
    )
