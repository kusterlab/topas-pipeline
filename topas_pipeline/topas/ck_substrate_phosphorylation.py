# compute substrate phosphorylation scores (a.k.a. kinase activity) for cytoplasmic
# kinases using the confident substrates from Florian Bayer's decryptM experiments

import sys
import logging
from pathlib import Path
from typing import Optional

import pandas as pd
import numpy as np
from tqdm import tqdm

import psite_annotation as pa

from .. import sample_metadata
from .. import utils
from . import expression_correction
from . import scoring
from ..preprocess import bridge_normalization
from ..preprocess import phospho_grouping

tqdm.pandas()

# hacky way to get the package logger instead of just __main__ when running as a module
logger = logging.getLogger(__package__ + "." + Path(__file__).stem)

VALIDATED_KINASES = [
    "JNK",
    "MTOR",
    "ERK",
    "CDK2",
    "ATM",
    "CDK9/12/13",
    "sMAP2K36",
    "sMAP2K12",
    "PLK1",
    "p38",
    "MAPKAPKs",
    "sCHEK1",
    "RSKs",
    "AKTs",
    "MSKs",
    "TBK1",
    "CAMK2s",
    "ROCK1/2",
    "CDK1",
    "GSK3",
    "S6Ks",
    "AURKB",
    "CDK4/6",
    "sAAK1",
    "SRCs",
    "ATR",
    "CHEK2",
    "FAKs",
    # "RIPK2",
    # "FRK",
    # "HERs",
    # "SGKs",
    # "PAK1/2/3",
    # "MAP3K20",
    # "CDK12/13",
    # "CDK9",
    # "sATR",
    # "sGSK3",
    # "sCDK1213",
    # "sCDK9",
    # "EPHAs",
]

META_COLS = [
    "Patient_Identifier",
    "Program",
    "code_oncotree",
    "breadcrumb_oncotree",
    "tissue_topology",
]


def calculate_cytoplasmic_kinase_scores(
    results_folder: str,
    sample_annotation_file: str,
    metadata_file: str,
    topas_kinase_substrate_file: str,
    topas_results_folder: str = "topas_scores",
    expression_corrected_input: bool = False,
    file_prefix = "ck_substrate_phosphorylation_",
    overwrite: bool = False,
):
    file_suffix = ""
    results_folder = Path(results_folder)
    metadata_file = Path(metadata_file)
    topas_kinase_substrate_file = Path(topas_kinase_substrate_file)

    if expression_corrected_input:
        expression_correction.correct_phospho_for_protein_expression(
            results_folder=results_folder, sample_annotation_file=sample_annotation_file
        )
        file_suffix = "_expressioncorrected"

    kinase_score_file = (
        results_folder
        / topas_results_folder
        / f"{file_prefix}scores{file_suffix}.tsv"
    )
    if kinase_score_file.is_file():
        if not overwrite:
            logger.info(
                f"Cytoplasmic kinase scoring skipped - found files already processed"
            )
            return
        logger.info(f"Found existing results but overwrite flag was set.")

    logger.info("Construct joint modified sequence groups between cohort and decryptM")
    df_patients = bridge_normalization.read_cohort_modified_sequence_groups(
        results_folder
    )
    automated_sites = load_decryptM_annotations(topas_kinase_substrate_file)
    peptidoform_groups = get_joint_modified_sequence_groups(
        df_patients, automated_sites
    )

    logger.info("Aggregate modified sequence groups in patient data")
    phospho = phospho_grouping.read_cohort_intensities_df(
        results_folder / f"preprocessed_pp2_agg_batchcorrected{file_suffix}.csv",
        sample_annotation_file=sample_annotation_file,
        keep_identification_metadata_columns=False,
    )
    pp_intensities_df = aggregate_patient_modified_sequence_groups(
        phospho, peptidoform_groups
    )

    logger.info("Aggregate modified sequence groups in decryptM data")
    automated_sites = aggregate_decryptm_modified_sequence_groups(
        automated_sites, peptidoform_groups
    )

    decryptM_kinases = get_patient_annotated_sites(pp_intensities_df, automated_sites)

    substrate_file = (
        results_folder
        / topas_results_folder
        / f"{file_prefix}peptide_intensities{file_suffix}.tsv"
    )
    write_substrate_peptides(
        pp_intensities_df,
        decryptM_kinases,
        substrate_file,
        kinases=VALIDATED_KINASES,
    )

    logger.info("Compute cytoplasmic kinase scores")
    centered_peptide_zvals_file = (
        results_folder
        / topas_results_folder
        / f"{file_prefix}peptides_centered{file_suffix}.tsv"
    )
    scores = compute_phosphorylation_scores(
        pp_intensities_df,
        decryptM_kinases,
        centered_peptide_zvals_file=centered_peptide_zvals_file,
        kinases=VALIDATED_KINASES,
    )

    save_scores(scores, kinase_score_file)
    if metadata_file is not None and len(str(metadata_file)) > 0:
        save_scores(scores, kinase_score_file, metadata_file)


def save_scores(
    scores: pd.DataFrame, kinase_score_file: Path, metadata_file: Optional[str] = None
):
    if metadata_file:
        metadata_df = sample_metadata.load(metadata_file)
        scores = merge_scores_with_sample_metadata(scores, metadata_df)
        kinase_score_file = kinase_score_file.with_name(
            kinase_score_file.stem + "_with_metadata.tsv"
        )

    scores.index.name = "Sample name"

    logger.info(f"Writing results to {kinase_score_file}")
    kinase_score_file.parent.mkdir(exist_ok=True)

    scores.to_csv(kinase_score_file, sep="\t", float_format="%.4f")


def load_scores(kinase_score_file: str) -> pd.DataFrame:
    return pd.read_csv(
        kinase_score_file,
        sep="\t",
        index_col=0,
    )


def merge_scores_with_sample_metadata(
    scores: pd.DataFrame, metadata_df: pd.DataFrame
) -> pd.DataFrame:
    """Merge scores DataFrame with metadata on sample name."""
    metadata_df["Sample name"] = utils.PATIENT_PREFIX + metadata_df["Sample name"]
    metadata_columns = metadata_df.columns.intersection(META_COLS).tolist()
    merged = scores.merge(
        right=metadata_df.set_index("Sample name")[metadata_columns],
        left_index=True,
        right_index=True,
        how="left",
    )

    score_cols = scores.columns.tolist()
    return merged[metadata_columns + score_cols]


def get_joint_modified_sequence_groups(
    df_patients: pd.DataFrame, df_decryptM: pd.DataFrame
) -> pd.DataFrame:
    df_annot = get_decryptm_modified_sequence_groups(df_decryptM)

    df_patients["Data set"] = "Patients"
    df_annot["Data set"] = "Annotated"
    df_combined = pd.concat([df_patients, df_annot])
    df_combined = explode_modified_sequence_groups(df_combined)

    # Aggregate duplicates
    df_combined = df_combined.groupby(["Modified sequence"])["Data set"].apply(set)
    df_combined = df_combined.reset_index()
    df_combined["dummy"] = 1  # dummy column for pa.aggregateModifiedSequenceGroups

    # Make the groups
    df_combined = pa.aggregateModifiedSequenceGroups(
        df_combined,
        experiment_cols=["dummy"],
        agg_cols={"Data set": lambda row: set.union(*row)},
        match_tolerance=2,
    )
    df_combined = df_combined[["Modified sequence group", "Data set"]]

    # Aggregate new Modified sequence group
    df_combined = df_combined.explode("Data set")

    peptidoform_groups = df_combined[["Modified sequence group"]]
    peptidoform_groups = peptidoform_groups.drop_duplicates()
    peptidoform_groups = explode_modified_sequence_groups(peptidoform_groups)
    return peptidoform_groups


def load_decryptM_annotations(
    topas_kinase_substrate_file: Path, filter_for_confident_relationships: bool = False
):
    df_decryptM = pd.read_csv(
        topas_kinase_substrate_file,
        sep="\t",
    )
    return df_decryptM


def get_decryptm_modified_sequence_groups(
    df_decryptM: pd.DataFrame,
) -> pd.DataFrame:
    df_decryptM = df_decryptM.rename(
        columns={
            "Modified sequence": "Modified sequence group",
        }
    )
    df_decryptM = (
        df_decryptM.groupby(["Modified sequence group"])
        .first()
        .reset_index()[["Modified sequence group"]]
    )
    return df_decryptM


def explode_modified_sequence_groups(df: pd.DataFrame) -> pd.DataFrame:
    df["Modified sequence"] = df["Modified sequence group"].str.split(";")
    df = df.explode("Modified sequence")
    return df


def load_phospho_data(
    results_folder: Path, pp_index: list[str], file_suffix: str = ""
) -> pd.DataFrame:
    phospho = pd.read_csv(
        results_folder / f"preprocessed_pp2_agg_batchcorrected{file_suffix}.csv",
        index_col=pp_index,
    )
    return phospho


def aggregate_patient_modified_sequence_groups(
    phospho: pd.DataFrame, peptidoform_groups: pd.DataFrame
) -> pd.DataFrame:
    phospho = 10**phospho
    pp_index = phospho.index.names
    phospho = phospho.reset_index()
    phospho = explode_modified_sequence_groups(phospho)

    phospho = (
        phospho.drop(columns="Modified sequence group")
        .merge(right=peptidoform_groups, on="Modified sequence", how="left")
        .drop(columns="Modified sequence")
    )

    # divide into numeric, string and groupby col and create agg dict
    numeric_cols = phospho.select_dtypes(include=np.number).columns
    group_col = "Modified sequence group"
    string_cols = [c for c in pp_index if c != group_col and c not in numeric_cols]
    agg_dict = {
        col: lambda x: ";".join(x.dropna().astype(str).unique()) for col in string_cols
    }
    for col in numeric_cols:
        agg_dict[col] = "mean"

    # Group by and aggregate --> then log10 transform
    aggregated_df = phospho.groupby("Modified sequence group", as_index=False).agg(
        agg_dict
    )
    aggregated_df[numeric_cols] = np.log10(aggregated_df[numeric_cols])
    aggregated_df = aggregated_df.set_index(pp_index)
    return aggregated_df


def aggregate_decryptm_modified_sequence_groups(
    automated_sites: pd.DataFrame, peptidoform_groups: pd.DataFrame
):
    # Regroup based on all peptidoforms
    automated_sites["Modified sequence"] = automated_sites[
        "Modified sequence"
    ].str.split(";")
    automated_sites = automated_sites.explode("Modified sequence")
    automated_sites = automated_sites.merge(
        right=peptidoform_groups, on="Modified sequence", how="left"
    ).drop(columns="Modified sequence")
    automated_sites = (
        automated_sites.groupby("Modified sequence group")["Kinase Families"]
        .apply(set)
        .apply(sorted)
        .str.join(";")
    )
    return automated_sites


def get_patient_annotated_sites(
    pp_agg_df: pd.DataFrame, automated_sites: pd.DataFrame
) -> pd.Series:
    patient_modified_sequence_groups = pp_agg_df.index.to_frame().reset_index(drop=True)
    decryptM_kinases = patient_modified_sequence_groups.merge(
        automated_sites, on="Modified sequence group", how="left"
    ).replace(np.nan, "")
    decryptM_kinases = decryptM_kinases.set_index(pp_agg_df.index.names)[
        "Kinase Families"
    ]
    return decryptM_kinases


def write_substrate_peptides(
    pp_intensities_df: pd.DataFrame,
    kinase_substrate_annotation_df: pd.DataFrame,
    substrate_file: Path,
    kinases: list[str] = None,
):
    exploded_substrates_df = explode_series(kinase_substrate_annotation_df)
    if not kinases:
        kinases = exploded_substrates_df.unique()

    substrate_modified_sequence_groups = (
        exploded_substrates_df[exploded_substrates_df.isin(kinases)]
        .index.get_level_values("Modified sequence group")
        .drop_duplicates()
    )

    substrate_intensities_df = (
        kinase_substrate_annotation_df.loc[substrate_modified_sequence_groups]
        .to_frame()
        .join(pp_intensities_df.loc[substrate_modified_sequence_groups], how="inner")
    )

    logger.info(f"Writing substrate intensities to {substrate_file}")
    substrate_file.parent.mkdir(exist_ok=True)
    substrate_intensities_df.to_csv(substrate_file, sep="\t", float_format="%.4g")


def compute_phosphorylation_scores(
    pp_intensities_df: pd.DataFrame,
    pp_annotation_series: pd.Series,
    centered_peptide_zvals_file: Path = None,
    patient_columns: pd.Index = None,
    kinases: list[str] = None,
    kinase_annot_level: str = "Kinase Families",
    explode: bool = True,
) -> pd.DataFrame:
    if explode:
        pp_annotation_series = explode_series(pp_annotation_series)

    if kinases:
        pp_annotation_series = pp_annotation_series[pp_annotation_series.isin(kinases)]

    if patient_columns is None:
        patient_columns = utils.filter_for_patient_columns(pp_intensities_df).columns

    scores: pd.DataFrame = pp_intensities_df.join(pp_annotation_series, how="inner")
    scores = scores.set_index(scores[kinase_annot_level], append=True).drop(
        columns=kinase_annot_level
    )

    # Scale input (dataframe) in mean or robust way
    center, _ = scoring.get_center_and_scale(
        scores[patient_columns],
        standardize=False,
        robust=False,
    )
    centered_vals = scores.sub(center, axis=0)

    if centered_peptide_zvals_file is not None:
        save_centered_peptide_zvals(
            centered_vals,
            centered_peptide_zvals_file,
            kinase_annot_level=kinase_annot_level,
        )

    centered_vals = centered_vals.reset_index()
    # drop pp index cols for aggregation
    centered_vals = centered_vals.drop(columns=pp_intensities_df.index.names)

    centered_vals = (
        centered_vals.groupby(by=pp_annotation_series.name)
        .progress_apply(
            scoring.z_aggregate,
            patient_columns=patient_columns,
            robust=False,
            standardize_output=True,
            agg_f="sum",
            clip_input=(-np.inf, np.inf),
            # clip_output=(-4, 4),
        )
        .T
    )
    return centered_vals


def determine_scoring_type(kinase_annot_level: str) -> str:
    if kinase_annot_level == "Kinase Families":
        return "ck_substrate_phosphorylation"
    if kinase_annot_level == "PSP Kinases":
        return "rtk_substrate_phosphorylation"
    else:
        return "protein_phosphorylation"


def save_centered_peptide_zvals(
    z_vals: pd.DataFrame, output_file: Path, kinase_annot_level: str = "Kinase Families"
):
    output_file.parent.mkdir(exist_ok=True)
    z_vals.to_csv(output_file, sep="\t", float_format="%.4f")


def explode_series(s: pd.Series, delimiter: str = ";") -> pd.Series:
    s: pd.Series = s.str.split(delimiter)
    return s.explode()


"""
python3 -m topas_pipeline.topas.ck_substrate_phosphorylation -c config_patients.json [--expression-corrected]
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
    parser.add_argument(
        "-e",
        "--expression-corrected",
        help="Use expression corrected phospho input.",
        action=argparse.BooleanOptionalAction,
    )
    args = parser.parse_args(sys.argv[1:])

    configs = config.load(args.config)

    calculate_cytoplasmic_kinase_scores(
        results_folder=configs.results_folder,
        sample_annotation_file=configs.sample_annotation,
        metadata_file=configs.metadata_annotation,
        topas_kinase_substrate_file=configs.clinic_proc.topas_kinase_substrate_file,
        expression_corrected_input=args.expression_corrected,
        overwrite=args.overwrite,
    )
