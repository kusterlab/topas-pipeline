import re
from typing import List
import logging

import numpy as np
import pandas as pd
from typing import List, Union

from .data_loader import DataLoader, extract_cohort_name

logger = logging.getLogger(__name__)


class LFQLoader(DataLoader):
    def __init__(self, evidence_files):
        self.evidence_files = evidence_files

    def load_data(self, use_cols: List[str]):
        """
        Convert LFQ evidence.txt to a pseudo TMT evidence.txt with a single TMT channel
        """
        # Make a small function so can be used with all loaders?
        all_batches = []
        for evidence_file_path in self.evidence_files:
            infile = open(evidence_file_path, "r")
            columns = infile.readline().split("\t")

            use_cols_tmp = [
                c for c in use_cols if not c.startswith("Reporter intensity corrected")
            ]

            if "Gene names" not in columns:
                use_cols_tmp.remove("Gene names")
            if "Fraction" not in columns:
                use_cols_tmp.remove("Fraction")

            # we used the data of prosit scored data where the picked group fdr was used together with prosit re-scoring
            evidence_df = pd.read_csv(
                evidence_file_path, usecols=use_cols_tmp, delimiter="\t"
            )

            evidence_df["Proteins"] = evidence_df["Proteins"].apply(
                keep_protein_name_from_fasta_header
            )
            evidence_df["Leading proteins"] = evidence_df["Leading proteins"].apply(
                keep_protein_name_from_fasta_header
            )

            # evidence_df = evidence_df.dropna(subset=['Proteins', 'Leading proteins'])

            # pretend that the LFQ data is TMT data with a single TMT channel
            evidence_df["Reporter intensity corrected 1"] = evidence_df["Intensity"]
            df = evidence_df

            # convert phospho modification notation, e.g. S(Phospho (STY)) => pS
            df["Modified sequence"] = df["Modified sequence"].str.replace(
                re.compile(r"([STY])\(Phospho \(STY\)\)"),
                lambda pat: f"p{pat.group(1)}",
                regex=True,
            )

            # we pretend each of the 247 experiments is its own batch
            cohort_name = extract_cohort_name(evidence_file_path)
            df["Batch"] = df["Experiment"].apply(lambda x: f"{cohort_name}_Batch{x}")
            all_batches.append(df)

        return all_batches


def keep_protein_name_from_fasta_header(value: Union[str, float]):
    if not pd.isnull(value):
        if "sp" in value:
            return value.split("|")[1]
        else:
            return value
    else:
        return np.nan
