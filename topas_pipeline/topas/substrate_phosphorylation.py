import os
import sys
import logging
import argparse
from pathlib import Path
from typing import Union

import pandas as pd

from .. import config
from . import scoring

logger = logging.getLogger(__name__)


def explode_psites_by_kinases(patient_df):
    """Filter for p-sites with a kinase annotation and explode by the kinase annotation.

    Args:
        patient_df: dataframe with a 'PSP Kinases' column that can contain multiple kinases.

    Returns:
        Exploded dataframe with a single kinase in the PSP Kinases column
    """
    psites_with_kinases = patient_df["PSP Kinases"].notna()
    patient_df.loc[psites_with_kinases, "PSP Kinases"] = patient_df.loc[
        psites_with_kinases, "PSP Kinases"
    ].str.split(";")
    patient_df = patient_df.loc[psites_with_kinases, :].explode("PSP Kinases")
    return patient_df


def add_kinase_to_kinase_list(kinase_list_sep: Union[str, float], kinase_to_add: str):
    if pd.isnull(kinase_list_sep):
        kinase_list = []
    else:
        kinase_list = kinase_list_sep.split(";")

    if kinase_to_add not in kinase_list:
        kinase_list.append(kinase_to_add)

    return ";".join(kinase_list)


def add_kinase_annotation(
    preprocessed_df: pd.DataFrame,
    modified_sequences: set[str],
    kinase_annotation: str,
):
    preprocessed_df.loc[
        preprocessed_df.index.isin(modified_sequences),
        "PSP Kinases",
    ] = preprocessed_df.loc[
        preprocessed_df.index.isin(modified_sequences),
        "PSP Kinases",
    ].apply(
        add_kinase_to_kinase_list, args=(kinase_annotation,)
    )
    return preprocessed_df


def compute_custom_kinase_scores(
    result_folder_path: str,
    kinase_results_output_folder: Union[str, Path],
    annotated_modified_sequences: dict[str, set[str]],
    extra_kinase_annot_bool: bool,
    force: bool = False,
) -> None:
    preprocessed_df = scoring.topas_score_preprocess(result_folder_path)
    preprocessed_df = preprocessed_df.drop(columns=["Unnamed: 0"])
    preprocessed_df = preprocessed_df.set_index("Modified sequence")

    for annotation_string, modified_sequences in annotated_modified_sequences.items():
        preprocessed_df = add_kinase_annotation(
            preprocessed_df, modified_sequences, annotation_string
        )

    preprocessed_df = preprocessed_df.reset_index()

    kinase_scoring(
        kinase_results_output_folder, preprocessed_df, extra_kinase_annot_bool, force
    )


def kinase_scoring(
    kinase_results_output_folder: Union[str, Path],
    preprocessed_df: pd.DataFrame,
    extra_kinase_annot_bool,
    force: bool = False,
):
    logger.info("Running kinase scoring module")

    if not force and os.path.exists(
        os.path.join(kinase_results_output_folder, "kinase_scores.tsv")
    ):
        logger.info(f"Kinase scoring skipped - found files already processed")
        return

    if not os.path.exists(kinase_results_output_folder):
        os.makedirs(kinase_results_output_folder)

    logger.info("  Calculate p-site weights")
    kinase_df = scoring.calculate_psite_weights(preprocessed_df)

    kinase_df = explode_psites_by_kinases(kinase_df)
    kinase_df.to_csv(
        os.path.join(
            kinase_results_output_folder, "kinase_annotated_patients_with_weights.tsv"
        ),
        sep="\t",
        index=False,
        float_format="%.4g",
    )

    logger.info("  Calculate modified sequence weights")
    kinase_summed_weights = scoring.calculate_modified_sequence_weights(
        kinase_df, "PSP Kinases"
    )

    # TODO: Export here; table with uncapped weights and uncapped zscores
    kinase_capped_values = scoring.cap_zscores_and_weights(kinase_summed_weights)

    logger.info("  Calculate weighted z-scores")
    kinase_scored_peptides = scoring.calculate_weighted_z_scores(kinase_capped_values)
    kinase_scored_peptides.to_csv(
        os.path.join(kinase_results_output_folder, "scored_peptides.tsv"),
        sep="\t",
        index=False,
        float_format="%.4g",
    )

    # TODO: get rid of this by using kinase names such as ERK_high_conf instead
    if extra_kinase_annot_bool:
        kinase_summed_weights_high_conf = scoring.calculate_modified_sequence_weights(
            kinase_df, "Kinase_high_conf"
        )
        kinase_capped_values_high_conf = scoring.cap_zscores_and_weights(
            kinase_summed_weights_high_conf
        )
        kinase_scored_peptides_high_conf = scoring.calculate_weighted_z_scores(
            kinase_capped_values_high_conf
        )
        kinase_scored_peptides_high_conf.to_csv(
            os.path.join(kinase_results_output_folder, "scored_peptides_high_conf.tsv"),
            sep="\t",
            index=False,
            float_format="%.4g",
        )

    logger.info("  Calculate kinase scores")
    kinase_first_level_scores = scoring.sum_weighted_z_scores(
        kinase_scored_peptides, by="PSP Kinases"
    )

    # scoring.plot_histograms_to_check_normality(kinase_first_level_scores)

    logger.info("  2nd level z-scoring, adding target space and writing results")
    kinase_scores = scoring.second_level_z_scoring(
        kinase_first_level_scores, "PSP Kinases"
    )
    kinase_spaces = scoring.get_target_space(
        annotated_peptides_df=kinase_df,
        scored_peptides_df=kinase_scored_peptides,
        grouping_by="PSP Kinases",
    )
    kinase_scores = pd.merge(
        left=kinase_spaces, right=kinase_scores, on="PSP Kinases", how="left"
    ).sort_values(by="PSP Kinases")
    kinase_scores = kinase_scores.set_index(["PSP Kinases", "No. of total targets"])
    kinase_scores.to_csv(
        os.path.join(kinase_results_output_folder, "kinase_scores.tsv"),
        sep="\t",
        float_format="%.4g",
    )

    if extra_kinase_annot_bool:
        kinase_first_level_scores_high_conf = scoring.sum_weighted_z_scores(
            kinase_scored_peptides_high_conf, by="Kinase_high_conf"
        )
        kinase_scores_high_conf = scoring.second_level_z_scoring(
            kinase_first_level_scores_high_conf, "Kinase_high_conf"
        )
        kinase_spaces_high_conf = scoring.get_target_space(
            annotated_peptides_df=kinase_df,
            scored_peptides_df=kinase_scored_peptides_high_conf,
            grouping_by="Kinase_high_conf",
        )
        kinase_scores_high_conf = pd.merge(
            left=kinase_spaces_high_conf,
            right=kinase_scores_high_conf,
            on="Kinase_high_conf",
            how="left",
        ).sort_values(by="Kinase_high_conf")
        kinase_scores_high_conf = kinase_scores_high_conf.set_index(
            ["Kinase_high_conf", "No. of total targets"]
        )
        kinase_scores_high_conf.to_csv(
            os.path.join(kinase_results_output_folder, "kinase_scores_high_conf.tsv"),
            sep="\t",
            float_format="%.4g",
        )
        # merge original kinase scores with high_conf and save
        kinase_scores = pd.concat([kinase_scores, kinase_scores_high_conf], axis=0)
        kinase_scores.index = kinase_scores.index.set_names(
            ["PSP Kinases", "No. of total targets"]
        )
        kinase_scores.to_csv(
            os.path.join(kinase_results_output_folder, "kinase_scores.tsv"),
            sep="\t",
            float_format="%.4g",
        )


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", help="Absolute path to configuration file.")
    parser.add_argument(
        "-k",
        "--kinase_results_folder",
        default="kinase_results",
        help="Relative path to kinase results folder inside the results folder.",
    )
    args = parser.parse_args(sys.argv[1:])
    configs = config.load(args.config)

    preprocessed_df = scoring.topas_score_preprocess(configs.results_folder)

    kinase_results_output_folder = os.path.join(
        configs.results_folder, args.kinase_results_folder
    )
    kinase_scoring(kinase_results_output_folder, preprocessed_df)
