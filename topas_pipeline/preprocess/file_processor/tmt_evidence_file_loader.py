import pandas as pd
from typing import List, Union

from topas_pipeline.data_loaders.data_loader import extract_batch_name
from topas_pipeline.io.reader import ReaderFactory
from topas_pipeline.preprocess.file_processor.evidence_file_loader import EvidenceFileLoader

class TMTEvidenceFileLoader(EvidenceFileLoader):
    def __init__(self, result_files: Union[str, List[str]], **kwargs):
        self.result_files = result_files
        self.options = kwargs

    def load_evidence_file(
        self, evidence_file_path: str, usecols: List[str] = None
    ) -> pd.DataFrame:
        reader = ReaderFactory.get_reader(evidence_file_path)
        df = reader.read(usecols=usecols)
        df = df.astype(
            {
                "Reverse": "category",
                "Experiment": "category",
                "Modifications": "category",
                "Potential contaminant": "category",
            }
        )
        df = self.modify_protein_and_peptide_info(df)
        df = self.filter_out_rows_with_no_reporter_ions(df)
        df = self.update_batch_info(df, evidence_file_path)
        return df
    
    @staticmethod
    def update_batch_info(df: pd.DataFrame, evidence_file_path: str) -> pd.DataFrame:
        batch_name = extract_batch_name(evidence_file_path)
        df["Batch"] = batch_name
        return df
