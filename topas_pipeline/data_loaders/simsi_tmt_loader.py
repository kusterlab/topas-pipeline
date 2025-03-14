import io
import re
from typing import List, Union
from pathlib import Path
import logging

import numpy as np
import pandas as pd

from .tmt_loader import TMTLoader
from .. import simsi
from .data_loader import (
    extract_batch_name,
    extract_experiment_name,
    test_batch_names_equals,
)

logger = logging.getLogger(__name__)


class SimsiTMTLoader(TMTLoader):

    def __init__(self, evidence_files, results_folder, simsi_folder, data_type):
        self.evidence_files = evidence_files
        self.results_folder = results_folder
        self.simsi_folder = simsi_folder
        self.data_type = data_type

    def load_data(self, use_cols: List[str]):
        simsi_evidence_file = simsi.find_simsi_evidence_file(
            self.results_folder, self.simsi_folder, self.data_type
        )
        all_batches = _parse_simsi_evidence_file(
            simsi_evidence_file, self.evidence_files, use_cols
        )
        return all_batches


def _parse_simsi_evidence_file(
    simsi_evidence_file: Union[str, Path, io.BytesIO],
    evidence_files: List[Union[str, Path]],
    use_cols: List[str],
) -> List[pd.DataFrame]:
    logger.info("Parsing SIMSI evidence file")

    # TODO: fix these columns in SIMSI?
    use_cols_new = list()
    for c in use_cols:
        if c == "Gene names":
            c = "Gene Names"
        elif c == "Potential contaminant":
            continue
        use_cols_new.append(c)

    use_cols_new.append("Transferred spectra count")
    df = pd.read_csv(simsi_evidence_file, sep="\t", usecols=use_cols_new)
    df.loc[df["Proteins"].str.contains("CON_", na=False), "Potential contaminant"] = "+"
    df = df.rename(columns={"Gene Names": "Gene names"})

    experiment_to_batch_name_dict = {
        extract_experiment_name(f): extract_batch_name(f) for f in evidence_files
    }
    df["Batch"] = df["Experiment"].replace(experiment_to_batch_name_dict)

    # Check that batches match in evidence file list and from simsi output
    test_batch_names_equals(experiment_to_batch_name_dict, df)

    # Change phosphosite notation in modified sequence
    df["Modified sequence"] = df["Modified sequence"].str.replace(
        re.compile(r"([STY])\(Phospho \(STY\)\)"),
        lambda pat: f"p{pat.group(1)}",
        regex=True,
    )

    df = df.replace(0, np.nan)
    df = df.loc[
        ~df.filter(regex="^Reporter intensity corrected").isnull().all(axis=1), :
    ]

    return [x for _, x in df.groupby("Batch")]
