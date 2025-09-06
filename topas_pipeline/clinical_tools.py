import logging
from itertools import compress
from pathlib import Path
from typing import List, Dict, Tuple, Union

import pandas as pd
import numpy as np
import psite_annotation as pa

from . import config
from .topas import annotation as topas_annotation

logger = logging.getLogger(__name__)

# # Old annotation file format column names
# TOPAS_SCORE_COLUMN = "basket"
# TOPAS_SCORE_COLUMNS = {
#     TOPAS_SCORE_COLUMN: "TOPAS annot",
# }
# TOPAS_SUBSCORE_COLUMN = "sub_basket"

TOPAS_SCORE_COLUMN = "TOPAS_score"
TOPAS_SCORE_COLUMNS = {
    TOPAS_SCORE_COLUMN: "TOPAS annot",
    "POI_category": "POI category",
    "POI_category_protein": "POI category (protein level)",
}
TOPAS_SUBSCORE_COLUMN = "TOPAS_subscore"
TOPAS_SUBSCORE_COLUMNS = {
    TOPAS_SUBSCORE_COLUMN: "TOPAS sublevel annot",
}

VALID_TOPAS_LEVELS = {
    "fp": ["expression", "kinase activity"],
    "pp": [
        "phosphorylation",
        "important phosphorylation",
        "kinase activity",
    ],
    "fp_with_ref": [
        "expression", "kinase activity"
    ],
    "pp_with_ref": [
        "phosphorylation", 
        "important phosphorylation",
        "kinase activity"
    ],
}


def add_phospho_annotations(
    df: pd.DataFrame,
    clinic_proc_config: config.ClinicProc,
) -> pd.DataFrame:
    """
    Phospho-site annotation of experimental data using in-house developed tool (MT) based mainly on Phosphosite Plus but also in vitro
    experiments from the lab of Ishihama.

    :param df: dataframe with measured peptide intensities
    :param clinic_proc_config: paths to files with PSP and TOPAS gene annotations
    """
    logger.info("Phosphosite annotation")

    df = df.reset_index()
    df = pa.addPeptideAndPsitePositions(
        df,
        clinic_proc_config.pspFastaFile,
        pspInput=True,
        returnAllPotentialSites=False,
    )

    df = pa.addPSPKinaseSubstrateAnnotations(
        df, clinic_proc_config.pspKinaseSubstrateFile, gene_name=True
    )
    df = pa.addPSPAnnotations(df, clinic_proc_config.pspAnnotationFile)
    df = pa.addPSPRegulatoryAnnotations(df, clinic_proc_config.pspRegulatoryFile)
    df["PSP_LT_LIT"] = df["PSP_LT_LIT"].apply(lambda x: max(x.split(";")))
    df["PSP_MS_LIT"] = df["PSP_MS_LIT"].apply(lambda x: max(x.split(";")))
    df["PSP_MS_CST"] = df["PSP_MS_CST"].apply(lambda x: max(x.split(";")))
    df.rename(
        columns={"Site positions": "Site positions identified (MQ)"}, inplace=True
    )
    df = pa.addPeptideAndPsitePositions(
        df, clinic_proc_config.pspFastaFile, pspInput=True, returnAllPotentialSites=True
    )
    df = df.set_index("Modified sequence", drop=True)

    # Add extra kinase annotations
    if len(str(clinic_proc_config.extra_kinase_annot)) > 0:
        df = add_extra_kinase_annotations(df, clinic_proc_config.extra_kinase_annot)

    return df


def get_extra_kinase_annot_df(path):
    extra_kinase_annot_df = pd.read_excel(path)
    extra_kinase_annot_df = extra_kinase_annot_df.rename(columns={'Gene name': 'Gene names', 'Site': 'Site positions', 'Modified_sequence': 'Modified sequence'})
    return extra_kinase_annot_df


def get_df_fill_in_missing_modified_sequence(df):
    df['Site positions'] = df['Site positions'].str.split(';')
    df = df.explode('Site positions').reset_index()
    df['Site positions'] = df['Site positions'].str.split('_').str[1]
    return df


def get_sites_to_get_peptides_for(extra_kinase_annot_df):
    # we need to get the peptides with these sites on
    sites_to_get_peptides_for = extra_kinase_annot_df.loc[
        extra_kinase_annot_df["Modified sequence"].isna(),
        ['Gene names', 'Site positions']
    ]
    return sites_to_get_peptides_for


def add_extra_kinase_annotations(
    df: pd.DataFrame, extra_kinase_annot: str
) -> pd.DataFrame:
    """
    Adds extra kinase annotations to the dataframe.

    :param df: dataframe with measured peptide intensities
    :param extra_kinase_annot: path to file with extra kinase annotations
    """
    logger.info("Extra kinase annotation")

    extra_kinase_annot_df = get_extra_kinase_annot_df(extra_kinase_annot)

    # get the sites we need to get peptides for
    sites_to_get_peptides_for = get_sites_to_get_peptides_for(extra_kinase_annot_df)

    # Prepare df by splitting 'Site positions' into lists
    df_exploded = get_df_fill_in_missing_modified_sequence(df)

    # subset to df to get the annot for (psite only given)
    df_to_fill = extra_kinase_annot_df[
        extra_kinase_annot_df['Modified sequence'].isna() &
        extra_kinase_annot_df['Gene names'].notna() &
        extra_kinase_annot_df['Site positions'].notna()
    ].copy()

    # loop over df_to_fill and fill in Modified sequence using the lookup
    for index, row in df_to_fill.iterrows():
        gene = row['Gene names']
        site = str(row['Site positions'])
        matches = df_exploded[
            (df_exploded['Gene names'].str.contains(gene)) &
            (df_exploded['Site positions'] == site)
        ]
        if not matches.empty:
            # If there are multiple matches, join them with ';'
            modified_sequences = ";".join(matches['Modified sequence'].unique())
            df_to_fill.at[index, 'Modified sequence'] = modified_sequences
            print(f"Filled Modified sequence for {gene} at site {site}: {modified_sequences}")
        else:
            print(f"No match found for {gene} at site {site}")

    df_to_fill = df_to_fill.assign(
        **{"Modified sequence": df_to_fill["Modified sequence"].str.split(";")}
    ).explode("Modified sequence")
    
    # and then we can merge back to extra_kinase_annot_df
    extra_kinase_annot_df = pd.concat([extra_kinase_annot_df[extra_kinase_annot_df['Modified sequence'].notna()], df_to_fill], ignore_index=True)

    # drop if any complete duplicate rows
    key_cols = ["Gene names", "Site positions", "Modified sequence"]
    extra_kinase_annot_df = (
        extra_kinase_annot_df
        .drop_duplicates(subset=key_cols)
        .reset_index(drop=True)
    )

    kinase_map = (
        extra_kinase_annot_df
        .groupby("Modified sequence")["PSP Kinase"]
        .apply(lambda x: ";".join(sorted(set(x))))
    )
    df["Kinase_high_conf"] = df.index.to_series().map(kinase_map)
    return df


def add_psp_urls(pp: pd.DataFrame) -> pd.DataFrame:
    """
    Function to add column to dataframe with URLs to proteins/isoforms that the
    peptides in the data belongs to:  https://www.phosphosite.org/. It uses already
    annotated information from PSP to check if any annotation exists. If it does, the URL
    is created from template and concatenated with the uniprot ID.

    :param pp: df to annotate with URL to PhosphoSitePlus
    :return: df with added annotation of URL to PhosphoSitePlus
    """
    pp[["PSP_URL", "PSP_URL_extra"]] = pp[["Start positions", "Proteins"]].apply(
        add_url_columns, axis=1, result_type="expand"
    )
    return pp


def add_url_columns(row) -> Tuple[str, List]:
    start_positions, proteins = row
    # create boolean list for (p-site, protein) pairs found in PSP or not
    # check for any modified peptide with start position different from -1
    found_psites = [
        int(value) > 0 for value in start_positions.split(";") if value != ""
    ]
    # row[0] is integer index of row and row[1] is column value
    proteins = list(compress(proteins.split(";"), found_psites))

    URLs = list()
    main_url = ""
    main_found = False

    # There can be found more than one main protein URL but then the first is used
    # and the rest goes with the isoform URLs
    url_start = "https://www.phosphosite.org/uniprotAccAction?id="
    for index, protein in enumerate(proteins):
        # don't allow isoforms (recognizable by "-" in their protein IDs) as main protein
        if "-" not in protein and not main_found:
            main_url = '=HYPERLINK("' + str(url_start) + str(protein) + '")'
            main_found = True
        else:
            URLs.append(str(url_start) + str(protein))

    return main_url, URLs


def add_topas_annotations(
    df: pd.DataFrame, annot_file: Union[str, Path], data_type: str, annot_type: str
) -> Tuple[pd.DataFrame, Dict]:
    """
    Adds columns with topas annotations and weights to a dataframe

    :param df: dataframe with a 'Gene names' column to be annotated
    :param annot_file: path to file with annotations
    :param data_type: either 'pp' for phospho or 'fp' for full proteome
    :param annot_type: either 'TOPAS_score', 'TOPAS_subscore', 'POI_category'
    """
    logger.info(f"Proteomics TOPAS annotation {data_type} {annot_type}")

    # Some dataframes might have empty cells so let's exchange them with nans
    df = df.replace(r"^\s*$", np.nan, regex=True)
    topas_annotation_df = topas_annotation.read_topas_annotations(annot_file)

    if "POI" in annot_type:
        if data_type == 'pp':
            
            if annot_type == 'POI_category_protein':
                gene_df = df[['Gene names']].fillna("")
                gene_df = gene_df['Gene names']
                annot_dict = create_identifier_to_topas_dict(topas_annotation_df, "POI", identifier_type="Gene names")
            else:
                gene_df = df.index.to_series()
                annot_dict = create_identifier_to_topas_dict(topas_annotation_df, "POI", identifier_type="Modified sequence")

            df[annot_type] = gene_df.apply(
                map_identifier_list_to_annot_types,
                annot_dict=annot_dict,
                annot_type=annot_type,
                with_weights=False,
            )
        else:
            if annot_type != 'POI_category_protein':
                gene_df = df.index.to_series()
                annot_dict = create_identifier_to_topas_dict(topas_annotation_df, "POI", identifier_type="Gene names")
                df[annot_type] = gene_df.apply(
                    map_identifier_list_to_annot_types,
                    annot_dict=annot_dict,
                    annot_type=annot_type,
                    with_weights=False,
                )

    if 'TOPAS' in annot_type:
        identifier_type = "Gene names"

        # Do we need to subset to expression because it's expression we annotate here? 
        topas_expression = topas_annotation_df[topas_annotation_df['level'] == 'expression']
        annot_dict = create_identifier_to_topas_dict(topas_annotation_df, data_type, identifier_type=identifier_type)
        
        if data_type == 'pp':
            gene_df = df[['Gene names']].fillna("") 

            # we give it another annot_type name to reach an POI if statement in the function
            # This makes us get the TOPAS_subscore category instead of TOPAS_score which is just 'other'
            df[annot_type] = gene_df['Gene names'].apply(
                map_identifier_list_to_annot_types,
                annot_dict=annot_dict,
                annot_type='TOPAS_POI',
                with_weights=False,
            )
        else:
            annot_dict = create_identifier_to_topas_dict(topas_expression, data_type, identifier_type=identifier_type)
            gene_df = df.index.to_series()        
            df[annot_type] = gene_df.apply(
                map_identifier_list_to_annot_types,
                annot_dict=annot_dict,
                annot_type=annot_type,
                with_weights=False,
            )

    return df


def map_identifier_list_to_annot_types(
    identifier_list: pd.Series,
    annot_dict: Dict[str, str],
    annot_type: str,
    with_weights: bool,
) -> pd.Series:
    """
    Takes a list of semicolon separated identifiers and returns the annot_levels
    matching the first identifier with annotations and weights
    Input identifier list has to be pd.Series and if method used via apply it has to be of type dataframe
    """
    # TODO: Create unit tests and refactor this function!
    annotations = []
    if isinstance(identifier_list, pd.Series):
        if ';' in identifier_list.iloc[0]:
            identifier_list = identifier_list.iloc[0].split(";")
    else:
        identifier_list = [identifier_list]

    for identifier in identifier_list:
        if identifier not in annot_dict:
            continue

        annot_type_in_column = TOPAS_SCORE_COLUMN
        if "POI" in annot_type:
            annot_type_in_column = TOPAS_SUBSCORE_COLUMN

        groups = annot_dict[identifier]["group"].split(";")
        annot_group = annot_dict[identifier][annot_type_in_column].split(";")
        annot_weight = annot_dict[identifier]["weight"].split(";")

        for i, group in enumerate(groups):
            annotations.append(annot_group[i])
    return ";".join(annotations)


def create_identifier_to_topas_dict(
    topas_annotation_df: pd.DataFrame,
    data_type: Union[str, None] = "fp",
    identifier_type: str = "Gene names",
) -> Dict[str, str]:
    """Collect all the TOPAS annotations per gene in a dictionary of {'gene_name': 'topas1;topas2;...'}

    Args:
        topas_annotation_df (pd.DataFrame): dataframe with TOPAS gene/p-peptide annotations
        data_type (Union[str, None], optional): One of "fp", "pp" and "POI". Defaults to "fp".
        identifier_type (str, optional): _description_. Defaults to "Gene names".

    Returns:
        Dict[str, str]: _description_
    """

    if "fp" in data_type:

        topas_annotation_df = topas_annotation_df[
            (topas_annotation_df["group"] != "OTHER")
            & (topas_annotation_df["level"].isin(VALID_TOPAS_LEVELS[data_type]))
        ]
    elif "pp" in data_type:
        # adjust above because we want to use == Other but subset to name TOPAS in the sub?
        topas_annotation_df = topas_annotation_df[
            (topas_annotation_df["group"] == "OTHER") 
            & (topas_annotation_df["TOPAS_subscore"].str.contains("paper"))  
        ]
    else:  # other proteins of interest (POI)
        topas_annotation_df = topas_annotation_df[
            topas_annotation_df["group"] == "OTHER"
        ]

    topas_annotation_df = topas_annotation_df.dropna(subset=[identifier_type])
    topas_annotation_df = topas_annotation_df.groupby([identifier_type]).agg(
        lambda x: ";".join(sorted(map(str, x)))
    )
    annot_dict = topas_annotation_df.to_dict("index")
    return annot_dict


"""
python3 -m topas_pipeline.clinical_tools -c config_patients.json -i <input_tsv> -o <output_tsv>
"""
if __name__ == "__main__":
    import argparse
    import time

    # Parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-c", "--config", required=True, help="Absolute path to configuration file."
    )
    parser.add_argument(
        "-i",
        "--input_file",
        required=True,
        help="Absolute path to a tab separated file.",
    )
    parser.add_argument(
        "-o", "--output_file", required=True, help="Absolute path to output file."
    )
    parser.add_argument(
        "-t",
        "--data_type",
        default="fp",
        help="Data type, either 'pp' or 'fp' (default: 'fp').",
    )

    args = parser.parse_args()

    configs = config.load(args.config)

    index_col = "Gene names"
    if args.data_type == "pp":
        index_col = "Modified sequence"

    df = pd.read_csv(args.input_file, sep="\t", index_col=index_col)
    annot_file = configs.clinic_proc.prot_baskets

    # Start pipeline
    t0 = time.time()
    for annot_type in ["TOPAS score", "POI"]:
        df = add_topas_annotations(
            df, annot_file, data_type=args.data_type, annot_type=annot_type
        )

    df.to_csv(args.output_file, sep="\t", float_format="%.6g")

    t1 = time.time()
    total = t1 - t0
    logger.info(f"TOPAS annotation finished in {total} seconds")
