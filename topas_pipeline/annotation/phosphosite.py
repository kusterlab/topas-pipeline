import logging
from itertools import compress
from pathlib import Path
from typing import List, Dict, Tuple, Union

import pandas as pd
import numpy as np
import psite_annotation as pa

from . import config

logger = logging.getLogger(__name__)


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
    return df
