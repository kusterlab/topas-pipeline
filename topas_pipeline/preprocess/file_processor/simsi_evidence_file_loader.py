import pandas as pd
from typing import List, Union

from topas_pipeline.io.reader import ReaderFactory
from topas_pipeline.preprocess.file_processor.evidence_file_loader import EvidenceFileLoader

class SimsiEvidenceFileLoader(EvidenceFileLoader):
    def __init__(self, result_files: Union[str, List[str]], **kwargs):
        self.result_files = result_files
        self.options = kwargs
    
    def load(self, usecols: List[str] = None) -> List[pd.DataFrame]:
        #Check if result_files is a string or list of len 1. If not raise ValueError
        if isinstance(self.result_files, list):
            if len(self.result_files) != 1:
                raise ValueError("SimsiEvidenceFileLoader only supports a single result file")
            else:
                evidence_file_path = self.result_files[0]
        else:
            evidence_file_path = self.result_files
        all_batches = self.load_evidence_file(evidence_file_path, usecols=usecols)
        return all_batches
        
    
    def load_evidence_file(
        self, evidence_file_path: str, usecols: List[str] = None
    ) -> pd.DataFrame:
        # TODO: fix these columns in SIMSI?
        usecols_new = list()
        for c in usecols:
            if c == "Gene names":
                c = "Gene Names"
            elif c == "Potential contaminant":
                continue
            usecols_new.append(c)

        usecols_new.append("Transferred spectra count")
        reader = ReaderFactory.get_reader(evidence_file_path)
        df = reader.read(usecols=usecols_new)
        df.loc[df["Proteins"].str.contains("CON_", na=False), "Potential contaminant"] = "+"
        df = df.rename(columns={"Gene Names": "Gene names"})
        df = self.modify_protein_and_peptide_info(df)
        df = self.filter_out_rows_with_no_reporter_ions(df)
        df = self.update_batch_info(df, self.options["sample_annotation_df"])
        # Check that batches match in evidence file list and from simsi output
        # test_batch_names_equals(experiment_to_batch_name_dict, df) # TODO: is this needed?
        return [x for _, x in df.groupby("Batch")]
    
    @staticmethod
    def update_batch_info(df: pd.DataFrame, sample_annotation_df: pd.DataFrame) -> pd.DataFrame:
        mapping = pd.DataFrame({
            "Batch_name": sample_annotation_df["Batch Name"].apply(lambda x: f"Batch{x}"),
            "Batch": sample_annotation_df[["Cohort","Batch Name"]].apply(lambda row: f"{row['Cohort']}_Batch{row['Batch Name']}", axis=1),
        }).set_index("Batch_name")["Batch"].to_dict()
        df["Batch"] = df["Experiment"].str.replace(r"(Batch\d+)_.*$", lambda m: mapping[m.group(1)], regex=True)
        return df
    