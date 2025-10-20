import os.path
import os
import sys
import numpy as np
import logging

import pandas as pd
from pathlib import Path
from typing import Union, Optional

from .. import metrics
from .. import utils
from .. import clinical_annotation
from .. import identification_metadata as id_meta
from . import annotation as topas_annotation

# hacky way to get the package logger instead of just __main__ when running as a module
logger = logging.getLogger(__package__ + "." + __file__)


def load_z_scores_fp(results_folder):
    z_scores_fp_df = metrics.read_measures(results_folder, "fp", [("z", "z-score")])[
        "z-score"
    ]
    z_scores_fp_df = z_scores_fp_df.rename(columns=lambda x: x.replace("zscore_", ""))
    z_scores_fp_df = utils.keep_only_sample_columns(z_scores_fp_df)

    # replace missing values but detected in batch with -4.0
    annot_fp = clinical_annotation.read_annotated_expression_file(results_folder, "fp")
    z_scores_fp_df = id_meta.replace_detected_in_batch(z_scores_fp_df, annot_fp, -4.0)

    return z_scores_fp_df


def load_z_scores_pp(results_folder):
    z_scores_pp_df = metrics.read_measures(results_folder, "pp", [("z", "z-score")])[
        "z-score"
    ]
    z_scores_pp_df = z_scores_pp_df.rename(columns=lambda x: x.replace("zscore_", ""))

    z_scores_pp_df = z_scores_pp_df.reset_index()
    z_scores_pp_df["Modified sequence"] = utils.short_phospho_notation(
        z_scores_pp_df["Modified sequence"]
    )
    z_scores_pp_df = z_scores_pp_df.set_index("Modified sequence")

    z_scores_pp_df = utils.keep_only_sample_columns(z_scores_pp_df)

    return z_scores_pp_df


@utils.validate_file_access
def load_protein_phosphorylation(
    results_folder,
    protein_results_folder: str = "protein_results",
    remove_multi_gene_groups: bool = False,
):
    protein_scores_file = os.path.join(
        results_folder, protein_results_folder, "protein_scores.tsv"
    )
    protein_phosphorylation_df = pd.read_csv(
        protein_scores_file,
        sep="\t",
        index_col="Gene names",
    )

    protein_phosphorylation_df = utils.keep_only_sample_columns(
        protein_phosphorylation_df
    )

    if remove_multi_gene_groups:
        # remove phosphoprotein groups which are a result of shared peptides
        protein_phosphorylation_df = protein_phosphorylation_df.loc[
            ~protein_phosphorylation_df.index.str.contains(";")
        ]
    return protein_phosphorylation_df


@utils.validate_file_access
def load_substrate_phosphorylation(
    results_folder,
    kinase_results_folder: str = "kinase_results",
    add_total_targets_to_index: bool = False,
):
    kinase_score_file = os.path.join(
        results_folder, kinase_results_folder, "kinase_scores.tsv"
    )
    index_col = ["PSP Kinases"]
    # TODO: this can probably be removed when refactoring report_creation.py
    if add_total_targets_to_index:
        index_col.append("No. of total targets")

    kinase_scores_df = pd.read_csv(
        kinase_score_file,
        sep="\t",
        index_col=index_col,
    )

    kinase_scores_df = utils.keep_only_sample_columns(kinase_scores_df)
    return kinase_scores_df


@utils.validate_file_access
def load_substrate_phosphorylation_num_targets(
    results_folder,
    kinase_results_folder: str = "kinase_results",
    add_total_targets_to_index: bool = False,
):
    kinase_score_file = os.path.join(
        results_folder, kinase_results_folder, "kinase_scores.tsv"
    )
    index_col = ["PSP Kinases"]
    if add_total_targets_to_index:
        index_col.append("No. of total targets")
    kinase_scores_df = pd.read_csv(
        kinase_score_file,
        sep="\t",
        index_col=index_col,
    )

    kinase_scores_df = kinase_scores_df.filter(regex=r"^targets_")
    return kinase_scores_df


def get_topas_subscore_calculator(
    z_scores_fp_df: pd.DataFrame,
    z_scores_pp_df: pd.DataFrame,
    protein_phosphorylation_df: pd.DataFrame,
    kinase_scores_df: pd.DataFrame,
):
    def calculate_topas_subscore(topas_subscore_annotation_df: pd.DataFrame):
        """Compute the TOPAS subscores using different rules and input data matrices.

        We cap the expression and protein phosphorylation between -4 and +4 to balance
        out the effect of gene amplifications. The substrate phosphorylation (=kinase
        score) is not capped because they are reliable indicators of aberrant RTK signaling.
        """
        scoring_rule = topas_subscore_annotation_df["Scoring rule"].iloc[0].lower()
        if scoring_rule == "highest z-score":
            topas_subscore = get_weighted_max(
                z_scores_fp_df,
                topas_subscore_annotation_df,
                "Gene names",
                "Gene names",
            )
        elif scoring_rule == "highest z-score (p-site)":
            topas_subscore = get_weighted_max(
                z_scores_pp_df,
                topas_subscore_annotation_df,
                "Modified sequence",
                "Modified sequence",
            )
        elif (
            scoring_rule
            == "highest protein phosphorylation score (2nd level z-score, fh)"
        ):
            topas_subscore = get_weighted_max(
                protein_phosphorylation_df,
                topas_subscore_annotation_df,
                "Gene names",
                "Gene names",
            )
        elif scoring_rule == "highest kinase score (2nd level z-score, fh)":
            topas_subscore = get_weighted_max(
                kinase_scores_df,
                topas_subscore_annotation_df,
                "PSP Kinases",
                "Gene names",
                clip_lower=None,
                clip_upper=None,
            )
        elif scoring_rule == "summed z-score":
            topas_subscore = get_summed_zscore(
                z_scores_fp_df,
                topas_subscore_annotation_df,
                "Gene names",
                "Gene names",
            )
        else:
            raise ValueError(f"Unknown scoring rule {scoring_rule}")

        return topas_subscore

    return calculate_topas_subscore


def get_weighted_max(
    z_score_df: pd.DataFrame,
    topas_subscore_df: pd.DataFrame,
    z_score_index_col: str,
    topas_subscore_index_col: str,
    clip_lower: Optional[int] = -4,
    clip_upper: Optional[int] = 4,
) -> pd.Series:

    z_scores = get_weighted_z_scores(
        z_score_df, topas_subscore_df, z_score_index_col, topas_subscore_index_col
    )
    z_scores = z_scores.clip(lower=clip_lower, upper=clip_upper)

    # take the maximum score per column (=sample)
    return z_scores.max()


def get_summed_zscore(
    z_score_df: pd.DataFrame,
    topas_subscore_df: pd.DataFrame,
    z_score_index_col: str,
    topas_subscore_index_col: str,
    clip_lower: Optional[int] = -4,
    clip_upper: Optional[int] = 4,
) -> pd.Series:

    z_scores = get_weighted_z_scores(
        z_score_df, topas_subscore_df, z_score_index_col, topas_subscore_index_col
    )
    z_scores = z_scores.clip(lower=clip_lower, upper=clip_upper)

    # calculate the summed score per column (=sample)
    return z_scores.sum()


def get_zscores(
    z_score_df: pd.DataFrame,
    topas_subscore_df: pd.DataFrame,
    z_score_index_col: str,
    topas_subscore_index_col: str,
):
    # only keep samples and z_score_index_col columns
    z_scores = z_score_df.astype(
        float
    )  # make sure all z-scores are floats and not objects
    z_scores = z_scores.reset_index()

    # deal with protein groups; in the z-score dataframe protein groups exists,
    # e.g. AKT1;AKT2;AKT3. We temporarily split them in separate rows to merge
    # with topas_subscore_df, where each row only contains a single protein. After
    # merging we combine those rows into a single row again.
    z_scores, z_score_index_col_exploded = utils.explode_on_separated_string(
        z_scores, z_score_index_col
    )
    # merge z-score dataframe with TOPAS genes
    z_scores = z_scores.merge(
        topas_subscore_df[[topas_subscore_index_col]],
        left_on=z_score_index_col_exploded,
        right_on=topas_subscore_index_col,
    )
    # merge proteins from the same protein group into a single row
    z_scores = z_scores.groupby(z_score_index_col).agg("first").reset_index()
    z_scores = z_scores.set_index(z_score_index_col)
    return z_scores


def get_weighted_z_scores(
    z_score_df: pd.DataFrame,
    topas_subscore_df: pd.DataFrame,
    z_score_index_col: str,
    topas_subscore_index_col: str,
) -> pd.DataFrame:
    """
    z_score_df is a pandas dataframe with z_scores as columns and genes/modified sequences as rows
    topas_subscore_df is a pandas dataframe with genes with weights as rows
    z_score_index_col is the column name with the identifier in z_score_df, e.g. "Gene names" or "Modified sequence"
    topas_subscore_index_col is the column name with the identifier in subbasket_df, e.g. "Gene names" or "Modified sequence"
    """
    # make sure all z-scores are floats and not objects
    z_scores = z_score_df.astype(float)
    z_scores = z_scores.reset_index()

    # deal with protein groups; in the z-score dataframe protein groups exists,
    # e.g. AKT1;AKT2;AKT3. We temporarily split them in separate rows to merge
    # with topas_subscore_df, where each row only contains a single protein. After
    # merging we combine those rows into a single row again.
    z_scores, z_score_index_col_exploded = utils.explode_on_separated_string(
        z_scores, z_score_index_col
    )

    # merge z-score dataframe with TOPAS genes
    z_scores = z_scores.merge(
        topas_subscore_df[[topas_subscore_index_col, "weight"]],
        left_on=z_score_index_col_exploded,
        right_on=topas_subscore_index_col,
        suffixes=(None, "_topas"),
    )
    if set(topas_subscore_df[topas_subscore_index_col]) != set(
        z_scores[z_score_index_col_exploded]
    ):
        logger.warning(
            f"Could not find all identifiers: {set(topas_subscore_df[topas_subscore_index_col]) - set(z_scores[z_score_index_col_exploded])}"
        )

    # if the index columns are the same a suffix is added in the merge
    if z_score_index_col == topas_subscore_index_col:
        topas_subscore_index_col += "_topas"

    # merge proteins from the same protein group into a single row
    z_scores = z_scores.groupby(z_score_index_col).agg("first").reset_index()

    z_scores = z_scores.set_index(z_score_index_col)
    # multiply the z-score by the associated weight (usually -1 or +1)
    weights = z_scores["weight"]
    z_scores = z_scores.drop(
        columns=[z_score_index_col_exploded, topas_subscore_index_col, "weight"]
    )
    z_scores = z_scores.multiply(weights, axis=0)
    return z_scores


def extract_topas_member_z_scores(
    results_folder: Union[str, Path],
    topas_annotation_file: Union[str, Path],
    topas_results_folder: str = "",
) -> None:
    """
    Extracts original p-site/gene z-scores / phosphorylation scores / kinase scores
    Requires that
    - [phospho|full_proteome]_z.csv files are generated by metrics.py
    - protein_results/protein_scores.tsv is generated by TOPAS_protein_phosphorylation_scoring.py
    - kinase_results/kinase_scores.tsv is generated by topas.substrate_phosphorylation.py
    """
    topas_annotation_df = topas_annotation.read_topas_annotations(topas_annotation_file)

    z_scores_fp_df = load_z_scores_fp(results_folder)
    z_scores_pp_df = load_z_scores_pp(results_folder)
    protein_phosphorylation_df = load_protein_phosphorylation(
        results_folder, remove_multi_gene_groups=True
    )
    kinase_scores_df = load_substrate_phosphorylation(results_folder)

    topas_annotation_df["weight"] = 1

    os.makedirs(
        os.path.join(results_folder, topas_results_folder, "basket_member_z_scores"),
        exist_ok=True,
    )

    for topas_name, topas_score_df in topas_annotation_df.groupby("TOPAS_score"):
        logger.info(f"Extracting TOPAS member z-scores for {topas_name}")
        z_score_dfs = []
        for topas_subscore_name, topas_subscore_df in topas_score_df.groupby(
            "TOPAS_subscore_level"
        ):
            scoring_rule = topas_subscore_df["Scoring rule"].iloc[0].lower()
            if scoring_rule == "highest z-score":
                topas_subscore_z_scores = get_weighted_z_scores(
                    z_scores_fp_df, topas_subscore_df, "Gene names", "Gene names"
                )
            elif scoring_rule == "highest z-score (p-site)":
                topas_subscore_z_scores = get_weighted_z_scores(
                    z_scores_pp_df,
                    topas_subscore_df,
                    "Modified sequence",
                    "Modified sequence",
                )
            elif (
                scoring_rule
                == "highest protein phosphorylation score (2nd level z-score, fh)"
            ):
                topas_subscore_z_scores = get_weighted_z_scores(
                    protein_phosphorylation_df,
                    topas_subscore_df,
                    "Gene names",
                    "Gene names",
                )
            elif scoring_rule == "highest kinase score (2nd level z-score, fh)":
                topas_subscore_z_scores = get_weighted_z_scores(
                    kinase_scores_df, topas_subscore_df, "PSP Kinases", "Gene names"
                )
            elif scoring_rule == "summed z-score":
                topas_subscore_z_scores = get_weighted_z_scores(
                    z_scores_fp_df, topas_subscore_df, "Gene names", "Gene names"
                )
            else:
                raise ValueError(f"Unknown scoring rule {scoring_rule}")

            topas_subscore_z_scores.index = topas_subscore_z_scores.index.map(
                lambda x: f"{x} - {topas_subscore_name}"
            )

            z_score_dfs.append(topas_subscore_z_scores)

        z_score_df = pd.concat(z_score_dfs).T
        z_score_df = z_score_df.drop(
            z_score_df[z_score_df.index.str.startswith("targets")].index
        )
        z_score_df.index = z_score_df.index.str.replace(r"-R\d{1}", "", regex=True)
        z_score_df.to_csv(
            os.path.join(
                results_folder,
                topas_results_folder,
                "basket_member_z_scores",
                f'basket_member_z_scores_{topas_name.replace("/", "_").replace(" ", "_")}.tsv',
            ),
            sep="\t",
            float_format="%.4g",
        )


if __name__ == "__main__":
    import argparse

    from .. import config

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-c", "--config", required=True, help="Absolute path to configuration file."
    )
    parser.add_argument(
        "-b",
        "--basket_results_folder",
        default="",
        help="Relative path to basket results folder inside the results folder.",
    )
    parser.add_argument(
        "-m",
        "--basket-member-zscores",
        action="store_true",
        help="For each basket, write original p-site/gene z-scores / phosphorylation scores / kinase scores to a file.",
    )
    args = parser.parse_args(sys.argv[1:])

    configs = config.load(args.config)
    # Create results folder and save configurations
    os.makedirs(
        os.path.join(configs.results_folder, args.basket_results_folder),
        exist_ok=True,
    )
    utils.init_file_logger(configs.results_folder, "TOPAS_scoring_log.txt")

    if args.basket_member_zscores:
        extract_topas_member_z_scores(
            configs.results_folder,
            topas_annotation_file=configs.clinic_proc.prot_baskets,
            topas_results_folder=args.basket_results_folder,
        )


def calculate_psite_weights(df: pd.DataFrame):
    """Calculates weight for each p-site according to the peptide occurrence.

    This function deals with the situation that a single p-site can be detected in
    multiple modified peptides. In the end, we want each p-site to only be weighted
    once, to prevent a single p-site from disproportionally contributing to the final
    TOPAS score.

    A peptide's occurrence is the number of patients the peptide was detected in. We
    trust the peptide (containing the p-site) with a higher occurrence more than those
    which were only detected sporadically.

    Args:
        df: DataFrame with z-scores in columns named 'pat_<patient_id>'

    Returns:
        Input dataframe with weight columns for each patient named 'weight_<patient_id>'
    """
    df = df.sort_values("Modified sequence")

    patient_zscore_columns = [col for col in df.columns if "pat_" in col]
    patient_weight_columns = [
        col.replace("pat_", "weight_") for col in patient_zscore_columns
    ]

    # calculate peptide count per patient, setting it to 0 where patient intensity is 0
    weights_df = df[patient_zscore_columns + ["Site positions"]]
    weights_df = weights_df.rename(columns=lambda x: x.replace("pat_", "weight_"))
    weights_df[patient_weight_columns] = (
        weights_df[patient_weight_columns].notna().mul(df["Peptide count"], axis=0)
    )

    # calculate psite count per patient
    psite_counts_df = weights_df.groupby("Site positions").sum()

    # divide peptide count by psite count per patient
    # TODO: check w unittest if the two steps also work
    # weights_df[patient_weight_columns] /= weights_df[['Site positions']].merge(psite_counts_df, on='Site positions', how='left', validate="many_to_one")[patient_weight_columns].values

    merged_df = weights_df[["Site positions"]].merge(
        psite_counts_df, on="Site positions", how="left", validate="many_to_one"
    )
    weights_df[patient_weight_columns] /= merged_df[patient_weight_columns].values

    weights_df = weights_df.drop(columns=["Site positions"])
    weights_df = weights_df.replace(0, np.nan)

    return pd.concat([df, weights_df], axis=1)


def calculate_modified_sequence_weights(
    patient_dataframe: pd.DataFrame, summing_column: str
):
    """Sums the p-site weights of each modified sequence constituting p-sites.

    Args:
        patient_dataframe: dataframe with p-site weights and z-scores
        summing_column: column name with the annotation that will eventually be grouped by

    Returns:
        dataframe with weights for each modified sequence
    """
    weight_dataframe = patient_dataframe.groupby([summing_column, "Modified sequence"])
    weight_dataframe = weight_dataframe.agg(
        **{
            weight_col: pd.NamedAgg(column=weight_col, aggfunc="sum")
            for weight_col in patient_dataframe.columns.tolist()
            if "weight_" in weight_col
        },
        **{
            pat_col: pd.NamedAgg(column=pat_col, aggfunc="first")
            for pat_col in patient_dataframe.columns.tolist()
            if "pat_" in pat_col
        },
    ).reset_index()
    return weight_dataframe


def cap_zscores_and_weights(patient_dataframe: pd.DataFrame) -> pd.DataFrame:
    """Cap all z-scores between -4 and 4 and all weights between 0 and 1.

    Args:
        patient_dataframe: dataframe with modified sequence weights and z-scores

    Returns:
        dataframe with capped modified sequence weights and z-scores
    """
    capped_dataframe = patient_dataframe.copy()
    capped_dataframe = cap_weights(capped_dataframe)
    capped_dataframe = cap_zscores(capped_dataframe)
    return capped_dataframe


def cap_weights(patient_dataframe: pd.DataFrame) -> pd.DataFrame:
    """Cap all all weights between 0 and 1.

    Args:
        patient_dataframe: dataframe with modified sequence weights

    Returns:
        dataframe with capped modified sequence weights
    """
    weightcols = [col for col in patient_dataframe.columns if "weight_" in col]
    patient_dataframe[weightcols] = (
        patient_dataframe[weightcols].astype(float).clip(upper=1)
    )
    return patient_dataframe


def cap_zscores(patient_dataframe: pd.DataFrame) -> pd.DataFrame:
    """Cap all z-scores between -4 and 4

    Args:
        patient_dataframe: dataframe with modified sequence z-scores

    Returns:
        dataframe with capped modified sequence z-scores
    """
    patcols = [col for col in patient_dataframe.columns if "pat_" in col]
    patient_dataframe[patcols] = (
        patient_dataframe[patcols].astype(float).clip(upper=4, lower=-4)
    )
    return patient_dataframe


def calculate_weighted_z_scores(df: pd.DataFrame):
    """Multiplies the patient z-scores by the p-site weight

    Args:
        df: DataFrame with z-scores in columns named 'pat_<patient_id>' and weights in columns named 'weight_<patient_id>'

    Returns:
        Input dataframe with added columns for each patient named 'weighted_<patient_id>'
    """
    patient_zscore_columns = [col for col in df.columns if "pat_" in col]
    patient_weight_columns = [
        col.replace("pat_", "weight_") for col in patient_zscore_columns
    ]

    weighted_zscores_df = df[patient_zscore_columns].mul(
        df[patient_weight_columns].rename(
            columns=lambda x: x.replace("weight_", "pat_")
        )
    )
    weighted_zscores_df = weighted_zscores_df.rename(
        columns=lambda x: x.replace("pat_", "weighted_")
    )

    return pd.concat([df, weighted_zscores_df], axis=1)


def sum_weighted_z_scores(patient_dataframe, by):
    relevant_columns = [
        col
        for col in patient_dataframe.columns
        if ("weight_" not in col) and ("pat_" not in col)
    ]
    score_dataframe = patient_dataframe[relevant_columns].groupby([by])
    # TODO: move aggregation into var for readability
    score_dataframe = score_dataframe.agg(
        **(
            {
                score: pd.NamedAgg(column=score, aggfunc="sum")
                for score in patient_dataframe.columns
                if "weighted_" in score
            }
        )
    ).reset_index()
    score_dataframe = score_dataframe.replace(0, np.nan)
    score_dataframe = score_dataframe.rename(
        columns={
            col: col.replace("weighted_", "pat_") for col in score_dataframe.columns
        }
    )
    return score_dataframe


def second_level_z_scoring(patient_dataframe, by_column, plot_histograms=False):
    # TODO: refactor this code
    patient_dataframe["mean"] = patient_dataframe[
        [
            col
            for col in patient_dataframe.columns
            if col not in [by_column, "Modified sequence"]
        ]
    ].mean(axis=1)
    patient_dataframe["stdev"] = patient_dataframe[
        [
            col
            for col in patient_dataframe.columns
            if col not in [by_column, "Modified sequence", "mean"]
        ]
    ].std(axis=1)
    for patient in patient_dataframe.columns:
        if patient in [by_column, "Modified sequence", "mean", "stdev"]:
            continue
        patient_dataframe[patient] = (
            patient_dataframe[patient] - patient_dataframe["mean"]
        ) / patient_dataframe["stdev"]
    return patient_dataframe


def calculate_peptide_occurrence(pp_df: pd.DataFrame):
    pp_df = pp_df.replace("nan", np.nan)
    pp_df = pp_df.replace("", np.nan)
    patient_list = pp_df.filter(regex="pat_").columns.tolist()
    pp_df = pp_df.astype({patient: "float32" for patient in patient_list})
    pp_df["Peptide count"] = pp_df[patient_list].notna().sum(axis=1)
    pp_df["Peptide occurrence"] = (
        pp_df["Peptide count"].astype(str) + "/" + str(len(patient_list))
    )
    return pp_df


def topas_score_preprocess(
    results_folder: str,
    concatenate_psp_and_extra_kinases: bool = False,
    discard_isoforms: bool = True,
):
    filepath = os.path.join(results_folder, "topas_score_preprocessed.tsv")

    if os.path.exists(filepath):
        logger.info(f"Reading previously generated file: {filepath}")
        return pd.read_csv(filepath, sep="\t")

    usecols = [
        "Modified sequence",
        "Gene names",
        "Site positions",
        "Site positions identified (MQ)",
        "PSP Kinases",
    ]
    if concatenate_psp_and_extra_kinases:
        usecols.append("Kinase_high_conf")

    annotations_df = pd.read_csv(
        os.path.join(results_folder, "annot_pp.csv"),
        index_col=0,
        usecols=usecols,
    )

    if concatenate_psp_and_extra_kinases:
        # combine comma separated strings of PSP Kinases and Kinase_high_conf
        annotations_df["PSP Kinases"] = (
            annotations_df["PSP Kinases"]
            .fillna("")
            .str.cat(annotations_df["Kinase_high_conf"].fillna(""), sep=";")
        )
        annotations_df["PSP Kinases"] = annotations_df["PSP Kinases"].str.strip(";")
        annotations_df["PSP Kinases"] = annotations_df["PSP Kinases"].replace(
            "", np.nan
        )

    annotations_df = annotations_df.rename(
        columns={
            "Site positions": "All site positions",
            "Site positions identified (MQ)": "Site positions",
        }
    )

    patients_zscores_df = pd.read_csv(
        os.path.join(results_folder, "phospho_measures_z.tsv"),
        keep_default_na=False,
        sep="\t",
    )
    patients_zscores_df.columns = patients_zscores_df.columns.str.removeprefix(
        "zscore_"
    )
    patient_columns = patients_zscores_df.filter(regex="pat_").columns.tolist()
    patients_zscores_df = patients_zscores_df.loc[
        :, ["Gene names", "Modified sequence"] + patient_columns
    ]

    patients_zscores_df = pd.merge(
        left=patients_zscores_df,
        right=annotations_df,
        on=["Modified sequence", "Gene names"],
        how="inner",
        validate="one_to_one",
    )
    patients_zscores_df = patients_zscores_df.set_index("Modified sequence")

    patients_zscores_df = calculate_peptide_occurrence(patients_zscores_df)

    patients_zscores_df = patients_zscores_df[
        patients_zscores_df["Site positions"].str.len() > 0
    ]
    patients_zscores_df["Site positions"] = patients_zscores_df[
        "Site positions"
    ].str.split(";")
    patients_zscores_df = patients_zscores_df.explode("Site positions")
    if discard_isoforms:
        # discards rows with isoform numbers in their Site position strings, e.g. P12345-2_S123
        patients_zscores_df = patients_zscores_df[
            ~patients_zscores_df["Site positions"].str.contains(
                r"^(?:[^\W_]+-\d+_[STY]\d+)$"
            )
        ]

    patients_zscores_df = patients_zscores_df.reset_index()
    patients_zscores_df.to_csv(filepath, sep="\t", float_format="%.4g")
    return patients_zscores_df


def calculate_per_patient_targets(scored_peptide_df: pd.DataFrame, grouping_by: str):
    patients = [i for i in scored_peptide_df.columns if "pat_" in i]
    per_patient_targets = (
        scored_peptide_df.groupby(grouping_by)[patients].count().reset_index()
    )
    per_patient_targets.rename(
        columns={
            col: col.replace("pat_", "targets_")
            for col in per_patient_targets.columns
            if "pat_" in col
        },
        inplace=True,
    )
    return per_patient_targets


def get_target_space(annotated_peptides_df, grouping_by, scored_peptides_df):
    total_targets = (
        annotated_peptides_df.groupby(grouping_by)["Modified sequence"]
        .nunique()
        .reset_index(name="No. of total targets")
    )
    per_patient_targets = calculate_per_patient_targets(scored_peptides_df, grouping_by)
    target_space = pd.merge(
        left=total_targets,
        right=per_patient_targets,
        on=grouping_by,
        how="inner",
        validate="1:1",
    )
    return target_space
