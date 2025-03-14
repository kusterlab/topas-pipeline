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
        use_cols = [
            c for c in use_cols if not c.startswith("Reporter intensity corrected")
        ]

        # Make a small function so can be used with all loaders?
        infile = open(self.evidence_files[0], "r")
        columns = infile.readline().split("\t")

        if "Gene names" not in columns:
            use_cols.remove("Gene names")
        if "Fraction" not in columns:
            use_cols.remove("Fraction")

        # we used the data of prosit scored data where the picked group fdr was used together with prosit re-scoring
        evidence_file = pd.read_csv(
            self.evidence_files[0], usecols=use_cols, delimiter="\t"
        )

        if "Gene names" not in columns and "|" in evidence_file["Proteins"][0]:
            evidence_file["Proteins"] = evidence_file["Proteins"].apply(
                keep_protein_name_from_fasta_header
            )
            evidence_file["Leading proteins"] = evidence_file["Leading proteins"].apply(
                keep_protein_name_from_fasta_header
            )

        # evidence_file = evidence_file.dropna(subset=['Proteins', 'Leading proteins'])

        # pretend that the LFQ data is TMT data with a single TMT channel
        evidence_file["Reporter intensity corrected 1"] = evidence_file["Intensity"]
        df = evidence_file.astype(
            {
                "Reverse": "category",
                "Experiment": "category",
                "Modifications": "category",
                "Potential contaminant": "category",
            }
        )

        # convert phospho modification notation, e.g. S(Phospho (STY)) => pS
        df["Modified sequence"] = df["Modified sequence"].str.replace(
            re.compile(r"([STY])\(Phospho \(STY\)\)"), lambda pat: f"p{pat.group(1)}"
        )

        # we pretend each of the 247 experiments is its own batch
        cohort_name = extract_cohort_name(self.evidence_files[0])
        df["Batch"] = df["Experiment"].apply(lambda x: f"{cohort_name}_Batch{x}")
        all_batches = [pd.DataFrame(v) for k, v in df.groupby("Batch")]

        return all_batches


def keep_protein_name_from_fasta_header(value: Union[str, float]):
    if not pd.isnull(value):
        if "sp" in value:
            return value.split("|")[1]
        else:
            return value
    else:
        return np.nan
