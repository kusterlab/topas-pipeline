# compute substrate phosphorylation scores (a.k.a. kinase activity) for receptor
# tyrosine kinases using the manually curated confident substrates by Annika Schneider and Firas Hamood

import sys
import logging
from pathlib import Path

import pandas as pd
import psite_annotation as pa

from .. import utils
from . import ck_substrate_phosphorylation
from ..preprocess import phospho_grouping

# hacky way to get the package logger instead of just __main__ when running as a module
logger = logging.getLogger(__package__ + "." + Path(__file__).stem)


def calculate_rtk_scores(
    results_folder: str,
    metadata_file: str,
    extra_kinase_annot: str,
    fasta_file: str,
    file_prefix = "rtk_substrate_phosphorylation_",
    topas_results_folder: str = "topas_scores",
    overwrite: bool = False,
):
    results_folder = Path(results_folder)
    kinase_score_file: Path = (
        results_folder / topas_results_folder / f"{file_prefix}scores.tsv"
    )
    if kinase_score_file.is_file():
        if not overwrite:
            logger.info(
                f"Receptor tyrosine substrate phosphorylation scoring skipped - found files already processed"
            )
            return
        logger.info(f"Found existing results but overwrite flag was set.")

    cohort_annotated_sites_df = get_annotated_modified_sequence_groups(
        results_folder / "preprocessed_pp.csv",
        extra_kinase_annot,
        fasta_file,
    )

    annotated_cohort_intensities_df = read_annotated_cohort_intensities_df(
        results_folder / "preprocessed_pp.csv",
        cohort_annotated_sites_df,
    )

    substrate_file = (
        results_folder / topas_results_folder / f"{file_prefix}peptide_intensities.tsv"
    )
    annotated_sites_mapping = get_annotated_sites_mapping(cohort_annotated_sites_df)
    ck_substrate_phosphorylation.write_substrate_peptides(
        annotated_cohort_intensities_df, annotated_sites_mapping, substrate_file
    )

    centered_peptide_zvals_file = (
        results_folder / topas_results_folder / f"{file_prefix}peptides_centered.tsv"
    )
    scores = ck_substrate_phosphorylation.compute_phosphorylation_scores(
        annotated_cohort_intensities_df,
        annotated_sites_mapping,
        centered_peptide_zvals_file=centered_peptide_zvals_file,
        kinase_annot_level="PSP Kinases",
        explode=True,
    )

    ck_substrate_phosphorylation.save_scores(scores, kinase_score_file)
    if metadata_file is not None and len(metadata_file) > 0:
        ck_substrate_phosphorylation.save_scores(
            scores, kinase_score_file, metadata_file
        )


def get_annotated_modified_sequence_groups(
    cohort_intensity_file: str,
    extra_kinase_annot: str,
    fasta_file: str,
) -> pd.DataFrame:
    logger.info("Reading cohort modified sequence groups")
    # loading only these three columns takes 20 sec instead of 4 minutes
    pp_df = pd.read_csv(
        cohort_intensity_file,
        usecols=["Gene names", "Modified sequence group", "Proteins"],
    )
    return annotate_modified_sequence_groups(pp_df, extra_kinase_annot, fasta_file)


def annotate_modified_sequence_groups(
    pp_df: pd.DataFrame,
    extra_kinase_annot: str,
    fasta_file: str,
) -> pd.DataFrame:
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


def explode_modified_sequence_group(df: pd.DataFrame) -> pd.DataFrame:
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
    cohort_intensity_file: str,
    cohort_annotated_sites_df: pd.DataFrame,
):
    logger.info("Read cohort intensities for annotated sites")
    # the column is called PSP Kinases because it was annotated with psite_annotation, it does contain the RTK substrates though
    keep_rows = cohort_annotated_sites_df.loc[
        cohort_annotated_sites_df["PSP Kinases"] != "", "index"
    ]
    skip_rows = cohort_annotated_sites_df.loc[
        ~cohort_annotated_sites_df["index"].isin(keep_rows), "index"
    ]
    skip_rows += 1  # add 1 to the row numbers to account for the header line

    return phospho_grouping.read_cohort_intensities_df(
        cohort_intensity_file, skiprows=skip_rows
    )


@utils.validate_file_access
def load_substrate_phosphorylation(
    results_folder,
    topas_results_folder: str = "topas_scores",
):
    kinase_score_file = (
        results_folder
        / topas_results_folder
        / "rtk_substrate_phosphorylation_scores.tsv"
    )

    kinase_scores_df = pd.read_csv(kinase_score_file, index_col=0, sep="\t")
    return prepare_kinase_scores_df(kinase_scores_df)


def prepare_kinase_scores_df(kinase_scores_df: pd.DataFrame) -> pd.DataFrame:
    kinase_scores_df = kinase_scores_df.T
    kinase_scores_df.index.name = "PSP Kinases"
    kinase_scores_df.index = kinase_scores_df.index.str.replace(
        "RTK-TOPAS", "TOPAS", regex=False
    )

    return kinase_scores_df


"""
python3 -m topas_pipeline.topas.rtk_substrate_phosphorylation -c config_patients.json
"""
if __name__ == "__main__":
    import argparse

    from .. import config

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-c", "--config", required=True, help="Absolute path to configuration file."
    )
    parser.add_argument(
        "-o",
        "--overwrite",
        action="store_true",
        help="Ignore existing results and recompute outputs.",
    )
    args = parser.parse_args(sys.argv[1:])

    configs = config.load(args.config)

    calculate_rtk_scores(
        results_folder=configs.results_folder,
        metadata_file=configs.metadata_annotation,
        extra_kinase_annot=configs.clinic_proc.extra_kinase_annot,
        fasta_file=configs.preprocessing.fasta_file,
        overwrite=args.overwrite,
    )
