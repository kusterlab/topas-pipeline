import logging
from pathlib import Path

import pandas as pd
import psite_annotation as pa

from .. import config
from ..topas import ck_substrate_phosphorylation as ck_scoring

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
    df = pa.addPeptideAndPsitePositions(
        df,
        clinic_proc_config.pspFastaFile,
        pspInput=True,
        returnAllPotentialSites=False,
    )

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
    df = pa.addPeptideAndPsitePositions(
        df, clinic_proc_config.pspFastaFile, pspInput=True, returnAllPotentialSites=True
    )
    df = df.set_index("Modified sequence", drop=True)
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
    df = df.merge(
        decryptM_annotations[["Modified sequence", "Kinase Families"]],
        on="Modified sequence",
        how="left",
    )
    df["Kinases (TOPAS)"] = (
        df[["Kinases (TOPAS)", "Kinase Families"]]
        .replace("", pd.NA)
        .agg(lambda x: ";".join(x.dropna()), axis=1)
    )
    df = df.drop(columns=["Kinase Families"])
    df = df.set_index("Modified sequence", drop=True)
    return df
