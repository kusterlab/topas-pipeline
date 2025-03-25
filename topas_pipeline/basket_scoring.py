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
from . import clinical_process
from . import identification_metadata as id_meta
from . import sample_metadata

# hacky way to get the package logger instead of just __main__ when running as a module
logger = logging.getLogger(__package__ + "." + __file__)


TUPAC_CATEGORIES = {
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


def get_tupac_rtks():
    return [tupac for tupac, category in TUPAC_CATEGORIES.items() if category == "RTK"]


def read_basket_scores(
    results_folder: Union[str, Path],
    data_type: str = "None",
    main_basket_only: bool = True,
    z_scored: bool = False,
) -> pd.DataFrame:
    """
    Read basket score results for report creation
    Requires one of [basket_scores_4th_gen.tsv|basket_scores.tsv]
    """
    z_scored_suffix = ""
    if z_scored:
        z_scored_suffix = "_zscored"

    for basket_score_file_name in [
        f"basket_scores_4th_gen{z_scored_suffix}.tsv",
        "basket_scores{z_scored_suffix}.tsv",
    ]:
        basket_scores_file_path = os.path.join(results_folder, basket_score_file_name)
        if not os.path.exists(basket_scores_file_path):
            continue

        try:
            basket_scores_df = pd.read_csv(
                basket_scores_file_path, sep="\t", index_col="Sample"
            )
        except PermissionError:
            raise PermissionError(
                f"Cannot open basket scores file, check if you have it open in Excel. {basket_scores_file_path}"
            )

        if basket_score_file_name == "basket_scores.tsv":
            if data_type != "None":
                try:
                    basket_scores_df = basket_scores_df.iloc[
                        :, basket_scores_df.columns.str.startswith(data_type.upper())
                    ]
                except ValueError:
                    raise ValueError(
                        f"Data type given for subset is not found in basket score file: {data_type}"
                    )
            # remove drug scores
            basket_scores_df = basket_scores_df.iloc[
                :, ~basket_scores_df.columns.str.contains("TOPAS")
            ]
            if main_basket_only:
                # remove RTK baskets
                basket_scores_df = basket_scores_df.iloc[
                    :, ~basket_scores_df.columns.str.contains("RTK")
                ]

        return basket_scores_df

    raise FileNotFoundError("No baskets score file found")


def read_sub_basket_scores(results_folder: Union[str, Path]) -> List:
    """
    Read basket score results for report creation
    Requires one of [basket_scores_4th_gen.tsv|basket_scores.tsv]
    """
    sub_basket_files = get_paths_to_sub_basket_files(results_folder)
    # read in df and combine
    list_dfs = []
    for file in sub_basket_files:
        sub_basket_file = os.path.join(results_folder, file)
        if os.path.exists(sub_basket_file):
            try:
                df = pd.read_csv(sub_basket_file, sep="\t", index_col="index")
            except PermissionError:
                raise PermissionError(
                    f"Cannot open subbasket scores file, check if you have it open in Excel. {sub_basket_file}"
                )
        else:
            raise FileNotFoundError(
                f"Subbaskets score file not found. {sub_basket_file}"
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
    subbasket_scores = pd.concat(list_dfs)
    return subbasket_scores


def get_paths_to_sub_basket_files(results_folder):
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
    annot_fp, _ = clinical_process.read_annotation_files(results_folder, False, "fp")
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
    baskets_file: Union[str, Path],
    data_types: List[str],
    basket_results_folder: str = "",
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
        os.path.join(results_folder, basket_results_folder, "basket_scores_4th_gen.tsv")
    ):
        logger.info(f"TOPAS scoring skipped - found files already preprocessed")
        return

    # TODO: find a way to read in only necessary columns - patients, topas annot+TOPAS sub annot --> not metadata
    metadata_df = sample_metadata.load(metadata_file)
    all_baskets_df = read_baskets_file_4th_gen(baskets_file)
    all_baskets_df = all_baskets_df[all_baskets_df["GROUP"] != "OTHER"]
    z_scores_fp_df = load_z_scores_fp(results_folder)
    z_scores_pp_df = load_z_scores_pp(results_folder)
    protein_phosphorylation_df = load_protein_phosphorylation(results_folder)
    kinase_scores_df = load_kinase_scores(results_folder)
    total_basket_scores = {}
    for basket_name, basket_df in all_baskets_df.groupby("TOPAS_SCORE"):
        logger.info(f"Calculating TOPAS scores for {basket_name}")
        subbasket_scores = {}
        for subbasket_name, subbasket_df in basket_df.groupby("TOPAS_subscore_level"):
            scoring_rule = subbasket_df["SCORING RULE"].iloc[0].lower()
            ligand = False
            if "Ligand" in subbasket_name:
                ligand = True
            if scoring_rule == "highest z-score":
                subbasket_score = get_weighted_max(
                    z_scores_fp_df, subbasket_df, "Gene names", "GENE NAME", ligand
                )
            elif scoring_rule == "highest z-score (p-site)":
                subbasket_score = get_weighted_max(
                    z_scores_pp_df,
                    subbasket_df,
                    "Modified sequence",
                    "MODIFIED SEQUENCE",
                    ligand,
                )
            elif (
                scoring_rule
                == "highest protein phosphorylation score (2nd level z-score, fh)"
            ):
                subbasket_score = get_weighted_max(
                    protein_phosphorylation_df,
                    subbasket_df,
                    "Gene names",
                    "GENE NAME",
                    ligand,
                )
            elif scoring_rule == "highest kinase score (2nd level z-score, fh)":
                subbasket_score = get_weighted_max(
                    kinase_scores_df, subbasket_df, "PSP Kinases", "GENE NAME", ligand
                )
            elif scoring_rule == "summed z-score":
                subbasket_score = get_summed_zscore(
                    z_scores_fp_df, subbasket_df, "Gene names", "GENE NAME", ligand
                )
            else:
                raise ValueError(f"Unknown scoring rule {scoring_rule}")
            subbasket_scores[f"{basket_name} - {subbasket_name}"] = subbasket_score

        subbasket_scores_df = pd.DataFrame.from_dict(subbasket_scores)

        # Use a mask. We only want to sum up positive subbasket scores to the total TUPAC score!
        subbasket_scores_df[f"{basket_name}_total_basket_score"] = (
            subbasket_scores_df.sum(axis=1)
        )
        total_basket_scores[basket_name] = subbasket_scores_df[
            f"{basket_name}_total_basket_score"
        ]

        # Remove replicate suffix before merging against the metadata table which only has sample names (i.e. without replicate suffixes)
        subbasket_scores_df["Sample name"] = subbasket_scores_df.index.str.replace(
            r"-R\d{1}", "", regex=True
        )
        if "Histologic subtype" in metadata_df.columns:
            subtype_df = metadata_df[["Sample name", "Histologic subtype"]]
        else:
            subtype_df = metadata_df[["Sample name"]]
        subbasket_scores_df = subbasket_scores_df.reset_index().merge(
            subtype_df, on="Sample name", how="left"
        )
        subbasket_output_file = os.path.join(
            results_folder,
            basket_results_folder,
            f'subbasket_scores_{basket_name.replace("/", "_").replace(" ", "_")}.tsv',
        )
        subbasket_scores_df.to_csv(subbasket_output_file, sep="\t", index=False, float_format="%.4g")

        # apply second-level z-scoring per basket (i.e. per column)
        # subbasket_scores_df = subbasket_scores_df.apply(zscore)
        # subbasket_output_file_zscored = os.path.join(results_folder, basket_results_folder, f'subbasket_scores_{basket_name.replace("/", "_").replace(" ", "_")}_zscored.tsv')
        # subbasket_scores_df.to_csv(subbasket_output_file_zscored, sep='\t', index=False)

        logger.info(
            f"Written basket results for {basket_name} to: {subbasket_output_file}"
        )

    basket_scores_df = pd.DataFrame.from_dict(total_basket_scores)
    save_basket_scores(
        basket_scores_df,
        os.path.join(
            results_folder, basket_results_folder, "basket_scores_4th_gen.tsv"
        ),
    )

    # apply second-level z-scoring per basket (i.e. per column)
    # basket_scores_df = basket_scores_df.apply(zscore)

    basket_scores_df = basket_scores_df.drop(
        basket_scores_df[basket_scores_df.index.str.startswith("targets")].index
    )
    measures = metrics.compute_measures(basket_scores_df.T)
    measures["z-score"].columns = measures["z-score"].columns.str.strip("zscore_")
    zscores = measures["z-score"].T

    save_rtk_scores_w_metadata(
        zscores, metadata_df, os.path.join(results_folder, "rtk_landscape.tsv")
    )
    save_basket_scores(
        zscores,
        os.path.join(
            results_folder, basket_results_folder, "basket_scores_4th_gen_zscored.tsv"
        ),
    )


def save_rtk_scores_w_metadata(
    basket_scores: pd.DataFrame, metadata_df: pd.DataFrame, out_file
):
    # subset to RTK using TUPAC_CATEGORIES keys and then map metadata

    basket_scores = basket_scores.loc[
        :, basket_scores.columns.isin(TUPAC_CATEGORIES.keys())
    ]
    if "code_oncotree" in metadata_df.columns:
        basket_scores = map_index_to_df(basket_scores, metadata_df)
        basket_scores = basket_scores.set_index("Sample name")

        # basket_scores = basket_scores.rename(index=lambda x: x.replace("score_", ""))
        basket_scores.index.name = "Sample name"
    basket_scores = basket_scores.fillna(0)
    basket_scores.to_csv(out_file, sep="\t", float_format="%.4g")


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
            .apply(count_basket_annotations)
        )

    # save ident and annot to file
    ident_annot = [item for item in dfs.keys() if "ident" in item or "annot" in item]
    ident_annot_df = [dfs[col] for col in ident_annot]

    ident_annot_df = pd.concat(ident_annot_df, axis=1)
    ident_annot_df.columns = ident_annot
    return ident_annot_df


def count_sample_identifications(df):
    return df.filter(regex="pat_", axis=1).mask(df == "").count()


def count_basket_annotations(column):
    temp = column.mask(column == "").dropna()
    return pd.Series(temp.index.where(temp.index != "", other=np.nan)).count()


def count_significant_baskets(
    df_zscored: pd.DataFrame, thresholds: List[float] = [2.0]
):
    """Computes number of significant basket scores per sample.

    Args:
        df_zscored: DataFrame with z-scores basket scores with one column per sample.
        thresholds: List of z-score thresholds to assess significance. Defaults to [2.0].

    Returns:
        Dataframe with sample names and number of significant basket scores.
    """
    significant_baskets_per_patient = (
        (df_zscored >= thresholds[0]).sum(axis=1).to_frame("num_significant_baskets")
    )
    if len(thresholds) > 1:
        for threshold in thresholds[1:]:
            significant_baskets_per_patient[
                f"num_significant_baskets_z>={threshold}"
            ] = (df_zscored >= threshold).sum(axis=1)
    return significant_baskets_per_patient


def basket_sheet_sanity_check(all_baskets_df: pd.DataFrame) -> None:
    """
    Explain

    TODO: add check that there are no protein groups in the gene names column, e.g. Gene1;Gene2
    """
    scoring_rules_found = set(all_baskets_df["SCORING RULE"].str.lower().unique())
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
        all_baskets_df[~all_baskets_df["MODIFIED SEQUENCE"].isnull()]["SCORING RULE"]
        .str.lower()
        .unique()
    )
    valid_scoring_rules = {"highest z-score (p-site)"}
    unknown_scoring_rules = list(scoring_rules_found - valid_scoring_rules)
    if len(unknown_scoring_rules) > 0:
        raise ValueError(
            f"Invalid scoring rule for entry with modified sequence: {unknown_scoring_rules}"
        )


def read_baskets_file_4th_gen(baskets_file):
    all_baskets_df = pd.read_excel(baskets_file)
    all_baskets_df = utils.whitespace_remover(all_baskets_df)
    for col in all_baskets_df.columns:
        if "WEIGHT" not in col:
            all_baskets_df[col] = all_baskets_df[col].str.strip()
    all_baskets_df["TOPAS_subscore_level"] = (
        all_baskets_df["TOPAS_SUBSCORE"] + " - " + all_baskets_df["LEVEL"]
    )
    all_baskets_df["WEIGHT"] = all_baskets_df["WEIGHT"].fillna(
        1
    )  # empty cell in WEIGHT column means weight = 1

    basket_sheet_sanity_check(all_baskets_df)

    return all_baskets_df


def get_weighted_max(
    z_score_df: pd.DataFrame,
    subbasket_df: pd.DataFrame,
    z_score_index_col: str,
    subbasket_index_col: str,
    ligand: bool,
) -> pd.Series:

    z_scores = get_weighted_z_scores(
        z_score_df, subbasket_df, z_score_index_col, subbasket_index_col
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
    subbasket_df: pd.DataFrame,
    z_score_index_col: str,
    subbasket_index_col: str,
    ligand: bool,
) -> pd.Series:

    z_scores = get_weighted_z_scores(
        z_score_df, subbasket_df, z_score_index_col, subbasket_index_col
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
    subbasket_df: pd.DataFrame,
    z_score_index_col: str,
    subbasket_index_col: str,
):
    # only keep samples and z_score_index_col columns
    z_scores = z_score_df.astype(
        float
    )  # make sure all z-scores are floats and not objects
    z_scores = z_scores.reset_index()

    # deal with protein groups; in the z-score dataframe protein groups exists,
    # e.g. AKT1;AKT2;AKT3. We temporarily split them in separate rows to merge
    # with subbasket_df, where each row only contains a single protein. After
    # merging we combine those rows into a single row again.
    z_scores, z_score_index_col_exploded = utils.explode_on_separated_string(
        z_scores, z_score_index_col
    )
    # merge z-score dataframe with basket genes
    z_scores = z_scores.merge(
        subbasket_df[[subbasket_index_col]],
        left_on=z_score_index_col_exploded,
        right_on=subbasket_index_col,
    )
    # merge proteins from the same protein group into a single row
    z_scores = z_scores.groupby(z_score_index_col).agg("first").reset_index()
    z_scores = z_scores.set_index(z_score_index_col)
    return z_scores


def get_weighted_z_scores(
    z_score_df: pd.DataFrame,
    subbasket_df: pd.DataFrame,
    z_score_index_col: str,
    subbasket_index_col: str,
) -> pd.DataFrame:
    """
    z_score_df is a pandas dataframe with z_scores as columns and genes/modified sequences as rows
    subbasket_df is a pandas dataframe with genes with weights as rows
    z_score_index_col is the column name with the identifier in z_score_df, e.g. "Gene names" or "Modified sequence"
    subbasket_index_col is the column name with the identifier in subbasket_df, e.g. "GENE NAME" or "Modified sequence"
    """
    # only keep samples and z_score_index_col columns
    z_scores = z_score_df.astype(
        float
    )  # make sure all z-scores are floats and not objects
    z_scores = z_scores.reset_index()

    # deal with protein groups; in the z-score dataframe protein groups exists,
    # e.g. AKT1;AKT2;AKT3. We temporarily split them in separate rows to merge
    # with subbasket_df, where each row only contains a single protein. After
    # merging we combine those rows into a single row again.
    z_scores, z_score_index_col_exploded = utils.explode_on_separated_string(
        z_scores, z_score_index_col
    )

    # merge z-score dataframe with basket genes
    z_scores = z_scores.merge(
        subbasket_df[[subbasket_index_col, "WEIGHT"]],
        left_on=z_score_index_col_exploded,
        right_on=subbasket_index_col,
    )
    if set(subbasket_df[subbasket_index_col]) != set(
        z_scores[z_score_index_col_exploded]
    ):
        logger.warning(
            f"Could not find all identifiers: {set(subbasket_df[subbasket_index_col]) - set(z_scores[z_score_index_col_exploded])}"
        )

    # merge proteins from the same protein group into a single row
    z_scores = z_scores.groupby(z_score_index_col).agg("first").reset_index()

    z_scores = z_scores.set_index(z_score_index_col)
    # multiply the z-score by the associated weight (usually -1 or +1)
    weights = z_scores["WEIGHT"]
    z_scores = z_scores.drop(
        columns=[z_score_index_col_exploded, subbasket_index_col, "WEIGHT"]
    )
    z_scores = z_scores.multiply(weights, axis=0)
    return z_scores


def calculate_basket_scores(
    measures: Dict[str, pd.DataFrame],
    basket_scores: pd.DataFrame,
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

        basket_scores = calculate_basket_scores_single(
            measures["score"], f"{data_type.upper()} - Scores", "basket", basket_scores
        )
        basket_scores = calculate_basket_scores_single(
            measures["score"], f"{data_type.upper()} - RTK Scores", "rtk", basket_scores
        )

    return basket_scores


def calculate_basket_scores_single(
    score_df: pd.DataFrame,
    sheet_name: str,
    basket_column: str,
    scores: pd.DataFrame,
    apply_z_score=False,
):
    logger.info(f"Calculating basket scores for {sheet_name}")

    patient_columns = [c for c in score_df.columns if c.startswith("score_")]

    # save number of identified annotated genes/p-sites for boxplot titles
    scores[f"{sheet_name}.num_identified"] = score_df[patient_columns].count()

    # filter only entries with a basket entry
    score_df = score_df[score_df[basket_column].notna()].copy(deep=False)

    # save number of identified annotated genes/p-sites for boxplot titles
    scores[f"{sheet_name}.num_annotated"] = score_df[patient_columns].count()

    score_df[basket_column] = score_df[basket_column].apply(lambda x: x.split(";"))
    score_df = score_df.astype({f"{basket_column}_weights": str}, errors="raise")
    score_df[f"{basket_column}_weights"] = score_df[f"{basket_column}_weights"].apply(
        lambda x: x.split(";")
    )

    # replicate each entry into multiple rows for each corresponding basket annotations
    score_df = score_df.explode([basket_column, f"{basket_column}_weights"])
    for basket, df_basket in score_df.groupby(basket_column):
        if "FP" in sheet_name:
            df_basket[patient_columns] = df_basket[patient_columns].multiply(
                df_basket[f"{basket_column}_weights"].astype(float), axis="index"
            )
        if (
            sheet_name == "FP - Scores"
            and basket_column == "basket"
            and basket in ["TK", "Immunotherapy"]
        ):
            scores[f"{sheet_name}.{basket}"] = df_basket[patient_columns].max(axis=0)
        else:
            scores[f"{sheet_name}.{basket}"] = df_basket[patient_columns].sum(axis=0)
        if apply_z_score:
            scores[f"{sheet_name}.{basket}"] = (
                scores[f"{sheet_name}.{basket}"]
                - scores[f"{sheet_name}.{basket}"].median()
            ) / scores[f"{sheet_name}.{basket}"].std()
    return scores


def save_basket_scores(basket_scores: pd.DataFrame, out_file: str):
    basket_scores = basket_scores.rename(index=lambda x: x.replace("score_", ""))
    basket_scores.index.name = "Sample"
    basket_scores = basket_scores.fillna(0)
    basket_scores.to_csv(out_file, sep="\t", float_format="%.4g")


def save_basket_scores_long_format(basket_scores, out_file):
    basket_names = basket_scores.columns[
        basket_scores.columns != "Sample"
    ].values.tolist()
    basket_scores_long = pd.melt(
        basket_scores.reset_index(), id_vars="Sample", value_vars=basket_names
    )
    basket_scores_long["Data type"] = [
        "FP" if "FP" in variable else "PP"
        for variable in basket_scores_long["variable"]
    ]
    basket_scores_long["Basket type"] = [
        "RTK" if "RTK Scores" in variable else "Main"
        for variable in basket_scores_long["variable"]
    ]
    basket_scores_long = basket_scores_long.loc[
        ~basket_scores_long["variable"].str.contains("TOPAS Drug"), :
    ]
    basket_scores_long["variable"] = (
        basket_scores_long["variable"].str.split(".").str.get(-1)
    )
    basket_scores_long.to_csv(out_file, sep="\t", float_format="%.4g")


def extract_basket_member_z_scores_4th_gen(
    results_folder: Union[str, Path],
    baskets_file: Union[str, Path],
    basket_results_folder: str = "",
) -> None:
    """
    Extracts original p-site/gene z-scores / phosphorylation scores / kinase scores
    Requires that
    - [phospho|full_proteome]_z.csv files are generated by metrics.py
    - protein_results/protein_scores.tsv is generated by TOPAS_protein_phosphorylation_scoring.py
    - kinase_results/kinase_scores.tsv is generated by TOPAS_kinase_scoring.py
    """

    all_baskets_df = read_baskets_file_4th_gen(baskets_file)

    z_scores_fp_df = load_z_scores_fp(results_folder)
    z_scores_pp_df = load_z_scores_pp(results_folder)
    protein_phosphorylation_df = load_protein_phosphorylation(results_folder)
    kinase_scores_df = load_kinase_scores(results_folder)

    all_baskets_df["WEIGHT"] = 1

    os.makedirs(
        os.path.join(results_folder, basket_results_folder, "basket_member_z_scores"),
        exist_ok=True,
    )

    for basket_name, basket_df in all_baskets_df.groupby("BASKET"):
        logger.info(f"Extracting basket member z-scores for {basket_name}")
        z_score_dfs = []
        for subbasket_name, subbasket_df in basket_df.groupby("subbasket_level"):
            scoring_rule = subbasket_df["SCORING RULE"].iloc[0].lower()
            if scoring_rule == "highest z-score":
                subbasket_z_scores = get_weighted_z_scores(
                    z_scores_fp_df, subbasket_df, "Gene names", "GENE NAME"
                )
            elif scoring_rule == "highest z-score (p-site)":
                subbasket_z_scores = get_weighted_z_scores(
                    z_scores_pp_df,
                    subbasket_df,
                    "Modified sequence",
                    "MODIFIED SEQUENCE",
                )
            elif (
                scoring_rule
                == "highest protein phosphorylation score (2nd level z-score, fh)"
            ):
                subbasket_z_scores = get_weighted_z_scores(
                    protein_phosphorylation_df, subbasket_df, "Gene names", "GENE NAME"
                )
            elif scoring_rule == "highest kinase score (2nd level z-score, fh)":
                subbasket_z_scores = get_weighted_z_scores(
                    kinase_scores_df, subbasket_df, "PSP Kinases", "GENE NAME"
                )
            elif scoring_rule == "summed z-score":
                subbasket_z_scores = get_weighted_z_scores(
                    z_scores_fp_df, subbasket_df, "Gene names", "GENE NAME"
                )
            else:
                raise ValueError(f"Unknown scoring rule {scoring_rule}")

            subbasket_z_scores.index = subbasket_z_scores.index.map(
                lambda x: f"{x} - {subbasket_name}"
            )

            z_score_dfs.append(subbasket_z_scores)

        z_score_df = pd.concat(z_score_dfs).T
        z_score_df = z_score_df.drop(
            z_score_df[z_score_df.index.str.startswith("targets")].index
        )
        z_score_df.index = z_score_df.index.str.replace(r"-R\d{1}", "", regex=True)
        z_score_df.to_csv(
            os.path.join(
                results_folder,
                basket_results_folder,
                "basket_member_z_scores",
                f'basket_member_z_scores_{basket_name.replace("/", "_").replace(" ", "_")}.tsv',
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
    utils.init_file_logger(configs.results_folder, "Basket_scoring_log.txt")

    compute_TOPAS_scores(
        configs.results_folder,
        configs.preprocessing.debug,
        data_types=configs.data_types,
        baskets_file=configs.clinic_proc.prot_baskets,
        metadata_file=configs.metadata_annotation,
        basket_results_folder=args.basket_results_folder,
    )

    if args.basket_member_zscores:
        extract_basket_member_z_scores_4th_gen(
            configs.results_folder,
            baskets_file=configs.clinic_proc.prot_baskets,
            basket_results_folder=args.basket_results_folder,
        )
