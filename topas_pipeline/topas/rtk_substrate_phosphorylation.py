# compute substrate phosphorylation scores (a.k.a. kinase activity) for receptor
# tyrosine kinases using the manually curated confident substrates by Annika Schneider and Firas Hamood

import sys
import logging
from pathlib import Path

import pandas as pd
import psite_annotation as pa

from topas_pipeline import sample_annotation
from topas_pipeline import preprocess_tools as prep
from . import ck_substrate_phosphorylation

# hacky way to get the package logger instead of just __main__ when running as a module
logger = logging.getLogger(__package__ + "." + Path(__file__).stem)


def calculate_rtk_scores(
    results_folder: str,
    metadata_file: str,
    extra_kinase_annot: str,
    sample_annotation_file: str,
    fasta_file: str,
):
    cohort_annotated_sites_df = get_annotated_modified_sequence_groups(
        results_folder,
        extra_kinase_annot,
        fasta_file,
    )

    annotated_sites_mapping = get_annotated_sites_mapping(cohort_annotated_sites_df)

    annotated_cohort_intensities_df = get_annotated_cohort_intensities_df(
        results_folder,
        sample_annotation_file,
        cohort_annotated_sites_df,
    )

    annotated_modified_sequence_groups = (
        annotated_cohort_intensities_df.index.get_level_values(
            "Modified sequence group"
        )
    )
    annotated_cohort_batch_corrected_df = get_annotated_cohort_batch_corrected_df(
        results_folder, annotated_modified_sequence_groups
    )

    annotated_cohort_intensities_df.update(annotated_cohort_batch_corrected_df)

    substrate_file = (
        f"{results_folder}/topas_scores/rtk_substrate_peptide_intensities.csv"
    )
    ck_substrate_phosphorylation.write_substrate_peptides(
        annotated_cohort_intensities_df, annotated_sites_mapping, substrate_file
    )

    scores = ck_substrate_phosphorylation.compute_substrate_phosphorylation_scores(
        annotated_cohort_intensities_df,
        annotated_sites_mapping,
    )

    kinase_score_file = (
        f"{results_folder}/topas_scores/rtk_substrate_phosphorylation_scores.csv"
    )
    ck_substrate_phosphorylation.save_scores_with_metadata_columns(
        scores, metadata_file, kinase_score_file
    )


def get_annotated_modified_sequence_groups(
    results_folder: str,
    extra_kinase_annot: str,
    fasta_file: str,
):
    logger.info("Reading cohort modified sequence groups")
    # loading only these three columns takes 20 sec instead of 4 minutes
    pp_df = pd.read_csv(
        f"{results_folder}/preprocessed_pp2_agg.csv",
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


def get_annotated_cohort_intensities_df(
    results_folder: str,
    sample_annotation_file: str,
    cohort_annotated_sites_df: pd.DataFrame,
):
    logger.info("Get cohort intensities for annotated sites")
    keep_rows = cohort_annotated_sites_df.loc[
        cohort_annotated_sites_df["PSP Kinases"] != "", "index"
    ]
    skip_rows = cohort_annotated_sites_df.loc[
        ~cohort_annotated_sites_df["index"].isin(keep_rows), "index"
    ]

    pp_df = pd.read_csv(
        f"{results_folder}/preprocessed_pp2_agg.csv", skiprows=skip_rows + 1
    )

    pp_df = pp_df.filter(regex=r"^(?!Identification metadata)")
    pp_df = pp_df.drop(
        columns=[
            "Modified sequence representative",
            "Modified sequence representative degree",
            "Delocalized sequence",
            "Proteins",
        ]
    )
    pp_df = pp_df.set_index(["Modified sequence group", "Gene names"])

    sample_annotation_df = pd.read_csv(sample_annotation_file)
    channel_to_sample_id_dict = sample_annotation.get_channel_to_sample_id_dict(
        sample_annotation_df,
        remove_qc_failed=True,
        remove_replicates=False,
    )

    index_cols = ["Modified sequence group", "Gene names"]
    pp_df = prep.rename_columns_with_sample_ids(
        pp_df.reset_index(),
        channel_to_sample_id_dict,
        index_cols=index_cols,
    )
    pp_df = pp_df.set_index(index_cols)
    return pp_df


def get_annotated_cohort_batch_corrected_df(
    results_folder: str, annotated_modified_sequence_groups: pd.Series
):
    logger.info("Get batch corrected cohort intensities for annotated sites")
    # only contains 6 out of ~50 annotated TOPAS-RTK sites
    df_patients = pd.read_csv(
        f"{results_folder}/preprocessed_pp2_agg_batchcorrected.csv",
        usecols=["Gene names", "Modified sequence group"],
    )
    df_patients = df_patients.reset_index()
    skip_rows_batch_corrected = df_patients.loc[
        ~df_patients["Modified sequence group"].isin(
            annotated_modified_sequence_groups
        ),
        "index",
    ]

    phospho_batch_corrected = pd.read_csv(
        f"{results_folder}/preprocessed_pp2_agg_batchcorrected.csv",
        index_col=[0, 1],
        skiprows=skip_rows_batch_corrected + 1,
    )
    return phospho_batch_corrected


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
