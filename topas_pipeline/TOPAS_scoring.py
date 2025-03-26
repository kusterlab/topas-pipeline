import os.path
import os
import sys
import numpy as np
import warnings
import logging

import pandas as pd
from pathlib import Path
from typing import List, Dict, Union

from . import metrics
from . import utils
from . import z_scoring as scoring
from . import clinical_annotation
from . import identification_metadata as id_meta
from . import sample_metadata

# hacky way to get the package logger instead of just __main__ when running as a module
logger = logging.getLogger(__package__ + "." + __file__)


TOPAS_CATEGORIES = {
    # the rules for the 5th generation
    "ALK": "RTK",
    "AXL": "RTK",
    "DDR1": "RTK",
    "DDR2": "RTK",
    "EGFR": "RTK",
    "ERBB2": "RTK",
    # "ERBB3": "RTK",
    "FGFR1": "RTK",
    "FGFR2": "RTK",
    "FGFR4": "RTK",
    "FLT1": "RTK",
    "FLT4": "RTK",
    # "INSR": "RTK",
    "KDR": "RTK",
    "KIT": "RTK",
    # "MERTK": "RTK",
    "NTRK2": "RTK",
    "NTRK3": "RTK",
    "PDGFRA": "RTK",
    "PDGFRB": "RTK",
    "IGF1R": "RTK",
    "EPHA2": "RTK",
    "RET": "RTK",
    "MET": "RTK",
    #'Proximal RTK signaling': 'downstream signaling',
}


def get_topas_rtks():
    return [topas_name for topas_name, category in TOPAS_CATEGORIES.items() if category == "RTK"]


def read_topas_scores(
    results_folder: Union[str, Path],
    data_type: str = "None",
    main_topas_only: bool = True,
    z_scored: bool = False,
) -> pd.DataFrame:
    """
    Read TOPAS score results for report creation
    Requires one of [basket_scores_4th_gen.tsv|basket_scores.tsv]
    """
    z_scored_suffix = ""
    if z_scored:
        z_scored_suffix = "_zscored"

    for topas_score_file_name in [
        f"basket_scores_4th_gen{z_scored_suffix}.tsv",
        "basket_scores{z_scored_suffix}.tsv",
    ]:
        topas_scores_file_path = os.path.join(results_folder, topas_score_file_name)
        if not os.path.exists(topas_scores_file_path):
            continue

        try:
            topas_scores_df = pd.read_csv(
                topas_scores_file_path, sep="\t", index_col="Sample"
            )
        except PermissionError:
            raise PermissionError(
                f"Cannot open TOPAS scores file, check if you have it open in Excel. {topas_scores_file_path}"
            )

        if topas_score_file_name == "basket_scores.tsv":
            if data_type != "None":
                try:
                    topas_scores_df = topas_scores_df.iloc[
                        :, topas_scores_df.columns.str.startswith(data_type.upper())
                    ]
                except ValueError:
                    raise ValueError(
                        f"Data type given for subset is not found in TOPAS score file: {data_type}"
                    )
            # remove drug scores
            topas_scores_df = topas_scores_df.iloc[
                :, ~topas_scores_df.columns.str.contains("TOPAS")
            ]
            if main_topas_only:
                # remove RTK TOPAS scores
                topas_scores_df = topas_scores_df.iloc[
                    :, ~topas_scores_df.columns.str.contains("RTK")
                ]

        return topas_scores_df

    raise FileNotFoundError("No TOPAS score file found")


def read_topas_subscores(results_folder: Union[str, Path]) -> pd.DataFrame:
    """
    Read TOPAS subscore results for report creation
    Requires one of [basket_scores_4th_gen.tsv|basket_scores.tsv]
    """
    topas_subscore_files = get_paths_to_topas_subscore_files(results_folder)
    # read in df and combine
    list_dfs = []
    for file in topas_subscore_files:
        topas_subscore_file = os.path.join(results_folder, file)
        if os.path.exists(topas_subscore_file):
            try:
                df = pd.read_csv(topas_subscore_file, sep="\t", index_col="index")
            except PermissionError:
                raise PermissionError(
                    f"Cannot open TOPAS subscores file, check if you have it open in Excel. {topas_subscore_file}"
                )
        else:
            raise FileNotFoundError(
                f"TOPAS subscore file not found. {topas_subscore_file}"
            )
        # get rid of total basket score column, sample name and sarcoma subtype
        # TODO: make this not dependent on lower and bigger case!
        df = df.drop(
            df.loc[
                :, df.columns.str.contains("total_basket_score|Sample|subtype|Subtype")
            ],
            axis=1,
        )
        list_dfs.append(df.transpose())
    list_dfs = [df.loc[:, ~df.columns.duplicated(keep="first")] for df in list_dfs]
    topas_subscores_df = pd.concat(list_dfs)
    return topas_subscores_df


def get_paths_to_topas_subscore_files(results_folder):
    files = [
        filename
        for filename in os.listdir(results_folder)
        if filename.startswith("subbasket_scores")
    ]
    return files


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


# TODO: merge with read_protein_scoring function in TOPAS_protein_phosphorylation_scoring.py
def load_protein_phosphorylation(results_folder):
    protein_phosphorylation_df = pd.read_csv(
        os.path.join(results_folder, "protein_results/protein_scores.tsv"),
        sep="\t",
        index_col="Gene names",
    )
    protein_phosphorylation_df = utils.keep_only_sample_columns(
        protein_phosphorylation_df
    )
    # remove phosphoprotein groups which are a result of shared peptides
    protein_phosphorylation_df = protein_phosphorylation_df.loc[
        ~protein_phosphorylation_df.index.str.contains(";")
    ]
    return protein_phosphorylation_df


# TODO: merge with read_kinase_scoring function in TOPAS_kinase_scoring.py
def load_kinase_scores(results_folder, kinase_results_folder: str = "kinase_results"):
    kinase_score_file = os.path.join(
        results_folder, kinase_results_folder, "kinase_scores.tsv"
    )

    kinase_scores_df = pd.read_csv(
        kinase_score_file,
        sep="\t",
        index_col="PSP Kinases",
    )
    kinase_scores_df = utils.keep_only_sample_columns(kinase_scores_df)
    return kinase_scores_df


def compute_TOPAS_scores(
    results_folder: Union[str, Path],
    debug: bool,
    metadata_file: Union[str, Path],
    topas_annotation_file: Union[str, Path],
    data_types: List[str],
    topas_results_folder: str = "",
) -> None:
    """
    Computes TOPAS subscores and TOPAS scores
    Requires that
    - [phospho|full_proteome]_z.csv files are generated by metrics.py
    - protein_results/protein_scores.tsv is generated by TOPAS_protein_phosphorylation_scoring.py
    - kinase_results/kinase_scores.tsv is generated by TOPAS_kinase_scoring.py
    """
    logger.info("Running TOPAS scoring module")
    if os.path.exists(
        os.path.join(results_folder, topas_results_folder, "basket_scores_4th_gen.tsv")
    ):
        logger.info(f"TOPAS scoring skipped - found files already preprocessed")
        return

    # TODO: find a way to read in only necessary columns - patients, topas annot+TOPAS sub annot --> not metadata
    metadata_df = sample_metadata.load(metadata_file)
    topas_annotation_df = read_topas_annotation_file(topas_annotation_file)
    topas_annotation_df = topas_annotation_df[topas_annotation_df["GROUP"] != "OTHER"]
    z_scores_fp_df = load_z_scores_fp(results_folder)
    z_scores_pp_df = load_z_scores_pp(results_folder)
    protein_phosphorylation_df = load_protein_phosphorylation(results_folder)
    kinase_scores_df = load_kinase_scores(results_folder)
    total_topas_scores = {}
    for topas_name, topas_score_df in topas_annotation_df.groupby("TOPAS_SCORE"):
        logger.info(f"Calculating TOPAS scores for {topas_name}")
        topas_subscores = {}
        for topas_subscore_name, topas_subscore_df in topas_score_df.groupby("TOPAS_subscore_level"):
            scoring_rule = topas_subscore_df["SCORING RULE"].iloc[0].lower()
            ligand = False
            if "Ligand" in topas_subscore_name:
                ligand = True
            if scoring_rule == "highest z-score":
                topas_subscore = get_weighted_max(
                    z_scores_fp_df, topas_subscore_df, "Gene names", "GENE NAME", ligand
                )
            elif scoring_rule == "highest z-score (p-site)":
                topas_subscore = get_weighted_max(
                    z_scores_pp_df,
                    topas_subscore_df,
                    "Modified sequence",
                    "MODIFIED SEQUENCE",
                    ligand,
                )
            elif (
                scoring_rule
                == "highest protein phosphorylation score (2nd level z-score, fh)"
            ):
                topas_subscore = get_weighted_max(
                    protein_phosphorylation_df,
                    topas_subscore_df,
                    "Gene names",
                    "GENE NAME",
                    ligand,
                )
            elif scoring_rule == "highest kinase score (2nd level z-score, fh)":
                topas_subscore = get_weighted_max(
                    kinase_scores_df, topas_subscore_df, "PSP Kinases", "GENE NAME", ligand
                )
            elif scoring_rule == "summed z-score":
                topas_subscore = get_summed_zscore(
                    z_scores_fp_df, topas_subscore_df, "Gene names", "GENE NAME", ligand
                )
            else:
                raise ValueError(f"Unknown scoring rule {scoring_rule}")
            topas_subscores[f"{topas_name} - {topas_subscore_name}"] = topas_subscore

        topas_subscores_df = pd.DataFrame.from_dict(topas_subscores)

        # Use a mask. We only want to sum up positive subbasket scores to the total TOPAS score!
        topas_subscores_df[f"{topas_name}_total_basket_score"] = (
            topas_subscores_df.sum(axis=1)
        )
        total_topas_scores[topas_name] = topas_subscores_df[
            f"{topas_name}_total_basket_score"
        ]

        # Remove replicate suffix before merging against the metadata table which only has sample names (i.e. without replicate suffixes)
        topas_subscores_df["Sample name"] = topas_subscores_df.index.str.replace(
            r"-R\d{1}", "", regex=True
        )
        if "Histologic subtype" in metadata_df.columns:
            subtype_df = metadata_df[["Sample name", "Histologic subtype"]]
        else:
            subtype_df = metadata_df[["Sample name"]]
        topas_subscores_df = topas_subscores_df.reset_index().merge(
            subtype_df, on="Sample name", how="left"
        )
        topas_subscores_output_file = os.path.join(
            results_folder,
            topas_results_folder,
            f'subbasket_scores_{topas_name.replace("/", "_").replace(" ", "_")}.tsv',
        )
        topas_subscores_df.to_csv(topas_subscores_output_file, sep="\t", index=False, float_format="%.4g")

        # apply second-level z-scoring per basket (i.e. per column)
        # subbasket_scores_df = subbasket_scores_df.apply(zscore)
        # subbasket_output_file_zscored = os.path.join(results_folder, basket_results_folder, f'subbasket_scores_{basket_name.replace("/", "_").replace(" ", "_")}_zscored.tsv')
        # subbasket_scores_df.to_csv(subbasket_output_file_zscored, sep='\t', index=False)

        logger.info(
            f"Written TOPAS results for {topas_name} to: {topas_subscores_output_file}"
        )

    topas_scores_df = pd.DataFrame.from_dict(total_topas_scores)
    save_topas_scores(
        topas_scores_df,
        os.path.join(
            results_folder, topas_results_folder, "basket_scores_4th_gen.tsv"
        ),
    )

    # apply second-level z-scoring per basket (i.e. per column)
    # basket_scores_df = basket_scores_df.apply(zscore)

    topas_scores_df = topas_scores_df.drop(
        topas_scores_df[topas_scores_df.index.str.startswith("targets")].index
    )
    measures = metrics.compute_measures(topas_scores_df.T)
    measures["z-score"].columns = measures["z-score"].columns.str.strip("zscore_")
    zscores = measures["z-score"].T

    save_rtk_scores_w_metadata(
        zscores, metadata_df, os.path.join(results_folder, "rtk_landscape.tsv")
    )
    save_topas_scores(
        zscores,
        os.path.join(
            results_folder, topas_results_folder, "basket_scores_4th_gen_zscored.tsv"
        ),
    )


def save_rtk_scores_w_metadata(
    topas_scores_df: pd.DataFrame, metadata_df: pd.DataFrame, out_file: str
) -> None:
    # subset to RTK using TOPAS_CATEGORIES keys and then map metadata

    topas_scores_df = topas_scores_df.loc[
        :, topas_scores_df.columns.isin(TOPAS_CATEGORIES.keys())
    ]
    if "code_oncotree" in metadata_df.columns:
        topas_scores_df = map_index_to_df(topas_scores_df, metadata_df)
        topas_scores_df = topas_scores_df.set_index("Sample name")

        # topas_scores_df = topas_scores_df.rename(index=lambda x: x.replace("score_", ""))
        topas_scores_df.index.name = "Sample name"
    topas_scores_df = topas_scores_df.fillna(0)
    topas_scores_df.to_csv(out_file, sep="\t", float_format="%.4g")


def map_index_to_df(df1, df2):
    # Step 1: Remove prefix "pre_" from the index of df1, if it exists
    df1_index = df1.index.str.replace("^pat_", "", regex=True)

    # Step 2: Use the cleaned index to merge df1 with df2 on 'sample' and 'batch'
    merged_df = df1.copy()
    merged_df.index = df1_index  # Replace index in df1 with cleaned index

    # Perform the mapping: merge df1 and df2 on the index and 'batch'
    result_df = merged_df.merge(
        df2[["Sample name", "code_oncotree"]],
        left_on=[merged_df.index],
        right_on=["Sample name"],
        how="left",
    )
    return result_df


def get_number_ident_annot_per_sample(
    results_folder: Union[str, Path], data_types: List[str]
) -> pd.DataFrame:
    """
    Calculate the number of identifications and annotations per sample based on the provided results folder and data types.

    Parameters:
    results_folder (Union[str, Path]): Path to the results folder.
    data_types (List[str]): List of data types to process.

    Returns:
    pd.DataFrame: DataFrame containing the number of identifications and annotations per sample.
    """
    dfs = dict()
    for data_type in data_types:
        if data_type == "fp":
            index_col = "Gene names"
            keep_default_na = True
        else:
            index_col = "Modified sequence"
            keep_default_na = False
        dfs[data_type] = pd.read_csv(
            os.path.join(results_folder, f"annot_{data_type}.csv"),
            index_col=index_col,
            keep_default_na=keep_default_na,
        )
        # atm PP is read with missing values instead of nan which is why masking is necessary. TODO: just fix this
        dfs[data_type + ".num_identified"] = count_sample_identifications(
            dfs[data_type]
        )
        tmt_channels = dfs[data_type].filter(regex="pat_").columns.tolist()
        dfs[data_type + ".num_annotated"] = (
            dfs[data_type][tmt_channels + ["basket"]]
            .set_index("basket")
            .apply(count_topas_annotations)
        )

    # save ident and annot to file
    ident_annot = [item for item in dfs.keys() if "ident" in item or "annot" in item]
    ident_annot_df = [dfs[col] for col in ident_annot]

    ident_annot_df = pd.concat(ident_annot_df, axis=1)
    ident_annot_df.columns = ident_annot
    return ident_annot_df


def count_sample_identifications(df):
    return df.filter(regex="pat_", axis=1).mask(df == "").count()


def count_topas_annotations(column):
    temp = column.mask(column == "").dropna()
    return pd.Series(temp.index.where(temp.index != "", other=np.nan)).count()


def count_significant_topas_scores(
    df_zscored: pd.DataFrame, thresholds: List[float] = [2.0]
):
    """Computes number of significant TOPAS scores per sample.

    Args:
        df_zscored: DataFrame with z-scores TOPAS scores with one column per sample.
        thresholds: List of z-score thresholds to assess significance. Defaults to [2.0].

    Returns:
        Dataframe with sample names and number of significant TOPAS scores.
    """
    significant_topas_per_patient = (
        (df_zscored >= thresholds[0]).sum(axis=1).to_frame("num_significant_baskets")
    )
    if len(thresholds) > 1:
        for threshold in thresholds[1:]:
            significant_topas_per_patient[
                f"num_significant_baskets_z>={threshold}"
            ] = (df_zscored >= threshold).sum(axis=1)
    return significant_topas_per_patient


def topas_sheet_sanity_check(topas_annotation_df: pd.DataFrame) -> None:
    """
    Explain

    TODO: add check that there are no protein groups in the gene names column, e.g. Gene1;Gene2
    """
    scoring_rules_found = set(topas_annotation_df["SCORING RULE"].str.lower().unique())
    valid_scoring_rules = {
        "highest z-score",
        "highest z-score (p-site)",
        "highest protein phosphorylation score (2nd level z-score, fh)",
        "highest kinase score (2nd level z-score, fh)",
        "summed z-score",
    }
    unknown_scoring_rules = list(scoring_rules_found - valid_scoring_rules)
    if len(unknown_scoring_rules) > 0:
        raise ValueError(f"Unknown scoring rules: {unknown_scoring_rules}")

    # validate that scoring rule is "highest z-score (p-site)" if modified sequence column is not empty
    scoring_rules_found = set(
        topas_annotation_df[~topas_annotation_df["MODIFIED SEQUENCE"].isnull()]["SCORING RULE"]
        .str.lower()
        .unique()
    )
    valid_scoring_rules = {"highest z-score (p-site)"}
    unknown_scoring_rules = list(scoring_rules_found - valid_scoring_rules)
    if len(unknown_scoring_rules) > 0:
        raise ValueError(
            f"Invalid scoring rule for entry with modified sequence: {unknown_scoring_rules}"
        )


def read_topas_annotation_file(topas_annotation_file: str):
    topas_annotation_df = pd.read_excel(topas_annotation_file)
    topas_annotation_df = utils.whitespace_remover(topas_annotation_df)
    for col in topas_annotation_df.columns:
        if "WEIGHT" not in col:
            topas_annotation_df[col] = topas_annotation_df[col].str.strip()
    topas_annotation_df["TOPAS_subscore_level"] = (
        topas_annotation_df["TOPAS_SUBSCORE"] + " - " + topas_annotation_df["LEVEL"]
    )
    topas_annotation_df["WEIGHT"] = topas_annotation_df["WEIGHT"].fillna(
        1
    )  # empty cell in WEIGHT column means weight = 1

    topas_sheet_sanity_check(topas_annotation_df)

    return topas_annotation_df


def get_weighted_max(
    z_score_df: pd.DataFrame,
    topas_subscore_df: pd.DataFrame,
    z_score_index_col: str,
    topas_subscore_index_col: str,
    ligand: bool,
) -> pd.Series:

    z_scores = get_weighted_z_scores(
        z_score_df, topas_subscore_df, z_score_index_col, topas_subscore_index_col
    )

    # cap the z-scores between -4 and +4 unless ligand
    if ligand:
        z_scores = z_scores.clip(lower=0, upper=4)
    else:
        z_scores = z_scores.clip(lower=-4, upper=4)

    # take the maximum score per column (=sample)
    return z_scores.max()


def get_summed_zscore(
    z_score_df: pd.DataFrame,
    topas_subscore_df: pd.DataFrame,
    z_score_index_col: str,
    topas_subscore_index_col: str,
    ligand: bool,
) -> pd.Series:

    z_scores = get_weighted_z_scores(
        z_score_df, topas_subscore_df, z_score_index_col, topas_subscore_index_col
    )
    # cap the z-scores between -4 and +4 unless ligand
    if ligand:
        z_scores = z_scores.clip(lower=0, upper=4)
    else:
        z_scores = z_scores.clip(lower=-4, upper=4)

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
    topas_subscore_index_col is the column name with the identifier in subbasket_df, e.g. "GENE NAME" or "Modified sequence"
    """
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
        topas_subscore_df[[topas_subscore_index_col, "WEIGHT"]],
        left_on=z_score_index_col_exploded,
        right_on=topas_subscore_index_col,
    )
    if set(topas_subscore_df[topas_subscore_index_col]) != set(
        z_scores[z_score_index_col_exploded]
    ):
        logger.warning(
            f"Could not find all identifiers: {set(topas_subscore_df[topas_subscore_index_col]) - set(z_scores[z_score_index_col_exploded])}"
        )

    # merge proteins from the same protein group into a single row
    z_scores = z_scores.groupby(z_score_index_col).agg("first").reset_index()

    z_scores = z_scores.set_index(z_score_index_col)
    # multiply the z-score by the associated weight (usually -1 or +1)
    weights = z_scores["WEIGHT"]
    z_scores = z_scores.drop(
        columns=[z_score_index_col_exploded, topas_subscore_index_col, "WEIGHT"]
    )
    z_scores = z_scores.multiply(weights, axis=0)
    return z_scores


def calculate_topas_scores(
    measures: Dict[str, pd.DataFrame],
    topas_scores_df: pd.DataFrame,
    data_type: str,
    score_method: str = "sum_weighted_zscore",
) -> pd.DataFrame:
    if score_method == "sum_weighted_zscore":
        measures["score"] = scoring.calculate_bounded_zscore_all(measures["z-score"])
    else:
        raise ValueError(f"score method '{score_method}' not yet supported")

    with warnings.catch_warnings():
        warnings.filterwarnings(
            action="ignore",
            category=pd.errors.PerformanceWarning,
            message=r".*DataFrame is highly fragmented*",
        )

        topas_scores_df = calculate_topas_scores_single(
            measures["score"], f"{data_type.upper()} - Scores", "basket", topas_scores_df
        )
        topas_scores_df = calculate_topas_scores_single(
            measures["score"], f"{data_type.upper()} - RTK Scores", "rtk", topas_scores_df
        )

    return topas_scores_df


def calculate_topas_scores_single(
    score_df: pd.DataFrame,
    sheet_name: str,
    topas_column: str,
    scores: pd.DataFrame,
    apply_z_score=False,
):
    logger.info(f"Calculating TOPAS scores for {sheet_name}")

    patient_columns = [c for c in score_df.columns if c.startswith("score_")]

    # save number of identified annotated genes/p-sites for boxplot titles
    scores[f"{sheet_name}.num_identified"] = score_df[patient_columns].count()

    # filter only entries with a TOPAS entry
    score_df = score_df[score_df[topas_column].notna()].copy(deep=False)

    # save number of identified annotated genes/p-sites for boxplot titles
    scores[f"{sheet_name}.num_annotated"] = score_df[patient_columns].count()

    score_df[topas_column] = score_df[topas_column].apply(lambda x: x.split(";"))
    score_df = score_df.astype({f"{topas_column}_weights": str}, errors="raise")
    score_df[f"{topas_column}_weights"] = score_df[f"{topas_column}_weights"].apply(
        lambda x: x.split(";")
    )

    # replicate each entry into multiple rows for each corresponding TOPAS annotations
    score_df = score_df.explode([topas_column, f"{topas_column}_weights"])
    for topas_name, topas_score_df in score_df.groupby(topas_column):
        if "FP" in sheet_name:
            topas_score_df[patient_columns] = topas_score_df[patient_columns].multiply(
                topas_score_df[f"{topas_column}_weights"].astype(float), axis="index"
            )
        if (
            sheet_name == "FP - Scores"
            and topas_column == "basket"
            and topas_name in ["TK", "Immunotherapy"]
        ):
            scores[f"{sheet_name}.{topas_name}"] = topas_score_df[patient_columns].max(axis=0)
        else:
            scores[f"{sheet_name}.{topas_name}"] = topas_score_df[patient_columns].sum(axis=0)
        if apply_z_score:
            scores[f"{sheet_name}.{topas_name}"] = (
                scores[f"{sheet_name}.{topas_name}"]
                - scores[f"{sheet_name}.{topas_name}"].median()
            ) / scores[f"{sheet_name}.{topas_name}"].std()
    return scores


def save_topas_scores(topas_scores_df: pd.DataFrame, out_file: str):
    topas_scores_df = topas_scores_df.rename(index=lambda x: x.replace("score_", ""))
    topas_scores_df.index.name = "Sample"
    topas_scores_df = topas_scores_df.fillna(0)
    topas_scores_df.to_csv(out_file, sep="\t", float_format="%.4g")


def save_topas_scores_long_format(topas_scores_df: pd.DataFrame, out_file: str):
    topas_names = topas_scores_df.columns[
        topas_scores_df.columns != "Sample"
    ].values.tolist()
    topas_scores_df_long = pd.melt(
        topas_scores_df.reset_index(), id_vars="Sample", value_vars=topas_names
    )
    topas_scores_df_long["Data type"] = [
        "FP" if "FP" in variable else "PP"
        for variable in topas_scores_df_long["variable"]
    ]
    topas_scores_df_long["Basket type"] = [
        "RTK" if "RTK Scores" in variable else "Main"
        for variable in topas_scores_df_long["variable"]
    ]
    topas_scores_df_long = topas_scores_df_long.loc[
        ~topas_scores_df_long["variable"].str.contains("TOPAS Drug"), :
    ]
    topas_scores_df_long["variable"] = (
        topas_scores_df_long["variable"].str.split(".").str.get(-1)
    )
    topas_scores_df_long.to_csv(out_file, sep="\t", float_format="%.4g")


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
    - kinase_results/kinase_scores.tsv is generated by TOPAS_kinase_scoring.py
    """

    topas_annotation_df = read_topas_annotation_file(topas_annotation_file)

    z_scores_fp_df = load_z_scores_fp(results_folder)
    z_scores_pp_df = load_z_scores_pp(results_folder)
    protein_phosphorylation_df = load_protein_phosphorylation(results_folder)
    kinase_scores_df = load_kinase_scores(results_folder)

    topas_annotation_df["WEIGHT"] = 1

    os.makedirs(
        os.path.join(results_folder, topas_results_folder, "basket_member_z_scores"),
        exist_ok=True,
    )

    for topas_name, topas_score_df in topas_annotation_df.groupby("BASKET"):
        logger.info(f"Extracting basket member z-scores for {topas_name}")
        z_score_dfs = []
        for topas_subscore_name, topas_subscore_df in topas_score_df.groupby("subbasket_level"):
            scoring_rule = topas_subscore_df["SCORING RULE"].iloc[0].lower()
            if scoring_rule == "highest z-score":
                topas_subscore_z_scores = get_weighted_z_scores(
                    z_scores_fp_df, topas_subscore_df, "Gene names", "GENE NAME"
                )
            elif scoring_rule == "highest z-score (p-site)":
                topas_subscore_z_scores = get_weighted_z_scores(
                    z_scores_pp_df,
                    topas_subscore_df,
                    "Modified sequence",
                    "MODIFIED SEQUENCE",
                )
            elif (
                scoring_rule
                == "highest protein phosphorylation score (2nd level z-score, fh)"
            ):
                topas_subscore_z_scores = get_weighted_z_scores(
                    protein_phosphorylation_df, topas_subscore_df, "Gene names", "GENE NAME"
                )
            elif scoring_rule == "highest kinase score (2nd level z-score, fh)":
                topas_subscore_z_scores = get_weighted_z_scores(
                    kinase_scores_df, topas_subscore_df, "PSP Kinases", "GENE NAME"
                )
            elif scoring_rule == "summed z-score":
                topas_subscore_z_scores = get_weighted_z_scores(
                    z_scores_fp_df, topas_subscore_df, "Gene names", "GENE NAME"
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

    from . import config

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

    compute_TOPAS_scores(
        configs.results_folder,
        configs.preprocessing.debug,
        data_types=configs.data_types,
        topas_annotation_file=configs.clinic_proc.prot_baskets,
        metadata_file=configs.metadata_annotation,
        topas_results_folder=args.basket_results_folder,
    )

    if args.basket_member_zscores:
        extract_topas_member_z_scores(
            configs.results_folder,
            topas_annotation_file=configs.clinic_proc.prot_baskets,
            topas_results_folder=args.basket_results_folder,
        )
