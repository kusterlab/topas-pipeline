import pandas as pd
from typing import List, Union

from topas_pipeline.data_loaders.data_loader import extract_cohort_name
from topas_pipeline.preprocess.file_processor.result_file_loader import (
    BaseResultFileLoader,
)
from topas_pipeline.io.reader import ReaderFactory


class EvidenceFileLoader(BaseResultFileLoader):
    def __init__(self, result_files: Union[str, List[str]], **kwargs):
        self.result_files = result_files
        self.options = kwargs

    def load(self, usecols: List[str] = None) -> List[pd.DataFrame]:
        if isinstance(self.result_files, str):
            self.result_files = [self.result_files]
        all_batches = []
        for evidence_file_path in self.result_files:
            df = self.load_evidence_file(evidence_file_path, usecols=usecols)
            all_batches.append(df)
        return all_batches

    def load_evidence_file(
        self, evidence_file_path: str, usecols: List[str] = None
    ) -> pd.DataFrame:
        reader = ReaderFactory.get_reader(evidence_file_path)
        df = reader.read(usecols=usecols)
        df = self.modify_protein_and_peptide_info(df)
        df = self.update_batch_info(df, evidence_file_path)
        return df

    @staticmethod
    def update_batch_info(df: pd.DataFrame, evidence_file_path: str) -> pd.DataFrame:
        # we pretend each of the 247 experiments is its own batch
        cohort_name = extract_cohort_name(evidence_file_path)
        df["Batch"] = df["Experiment"].apply(lambda x: f"{cohort_name}_Batch{x}")
        return df
