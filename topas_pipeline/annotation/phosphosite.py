import logging
from pathlib import Path

import pandas as pd
import psite_annotation as pa

from .. import config
from .. import utils
from ..topas import ck_substrate_phosphorylation as ck_scoring
from ..preprocess import phospho_grouping

logger = logging.getLogger(__name__)


def add_phospho_annotations(
    df: pd.DataFrame,
    clinic_proc_config: config.ClinicProc,
) -> pd.DataFrame:
    """
    Phospho-site annotation of experimental data using in-house developed tool (MT) based mainly on Phosphosite Plus

    :param df: dataframe with measured peptide intensities
    :param clinic_proc_config: paths to files with PSP and TOPAS gene annotations
    """
    logger.info("Phosphosite annotation")

    df = df.reset_index()

    # add all potential phosphosites on each sequence
    df = pa.addPeptideAndPsitePositions(
        df, clinic_proc_config.pspFastaFile, pspInput=True, returnAllPotentialSites=True
    )
    df.rename(columns={"Site positions": "Site positions (PSP)"}, inplace=True)

    df = pa.addPeptideAndPsitePositions(
        df,
        clinic_proc_config.pspFastaFile,
        pspInput=True,
        returnAllPotentialSites=False,
    )

    # add semicolon to columns which will be concatenated. this allows us to use the
    # fast "sum" aggfunc instead of a slow custom string function
    concat_cols = ["Site positions", "Site sequence context", "Kinase Families"]
    df[concat_cols] = df[concat_cols].apply(lambda column: column.astype(str) + ";")

    # the PSP annotation functions can deal with semicolon separated strings 
    # in the "Site positions" column, so we groupby the modified sequence groups
    # before annotating.
    df = (
        df.groupby(by=phospho_grouping.INDEX_COLS)
        .agg(
            {
                "Matched proteins": "first",
                "Start positions": "first",
                "End positions": "first",
                "Site positions (PSP)": "first",
                "Site positions": "sum",
                "Site sequence context": "sum",
                "Kinase Families": "sum",
            }
        )
        .reset_index()
    )

    for col in concat_cols:
        df[col] = df[col].apply(utils.csv_unique)

    df = pa.addPSPKinaseSubstrateAnnotations(df, clinic_proc_config.extra_kinase_annot)
    df = df.rename(columns={"PSP Kinases": "Kinases (TOPAS)"})

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
    df = df.set_index("Modified sequence group", drop=True)
    return df


def add_ck_substrate_annotations(
    df: pd.DataFrame, results_folder: str, topas_kinase_substrate_file: str
):
    joint_peptidoform_groups = ck_scoring.get_joint_modified_sequence_groups(
        Path(results_folder), topas_kinase_substrate_file
    )

    decryptm_peptidoforms = ck_scoring.load_decryptM_annotations(
        topas_kinase_substrate_file
    )
    decryptm_peptidoform_groups = (
        ck_scoring.aggregate_decryptm_modified_sequence_groups(
            decryptm_peptidoforms, joint_peptidoform_groups
        )
    )

    decryptM_annotations = joint_peptidoform_groups.merge(
        decryptm_peptidoform_groups, on="Modified sequence group"
    )
    df = df.reset_index().merge(
        decryptM_annotations[["Modified sequence", "Kinase Families"]],
        on="Modified sequence",
        how="left",
    )
    df["Kinase Families"] = df["Kinase Families"].fillna("")
    df = df.set_index("Modified sequence", drop=True)
    return df


def combine_topas_annotations(df: pd.DataFrame) -> pd.DataFrame:
    df["Kinases (TOPAS)"] = (
        df[["Kinase Families", "Kinases (TOPAS)"]]
        .replace("", pd.NA)
        .agg(lambda x: ";".join(x.dropna()), axis=1)
    )
    df = df.drop(columns=["Kinase Families"])
    return df
