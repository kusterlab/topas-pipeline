from abc import ABC, abstractmethod
import pandas as pd

import topas_pipeline.preprocess.file_processor.utils as file_processor_utils


class BaseResultFileLoader(ABC):
    @abstractmethod
    def load(self):
        raise NotImplementedError(
            f"Class {self.__class__.__name__} has not implemented the load() method"
        )

    @staticmethod
    def modify_protein_and_peptide_info(df: pd.DataFrame) -> pd.DataFrame:
        return file_processor_utils.modify_protein_and_peptide_info(df)
    
    @staticmethod
    def add_reporter_intensity_corrected_1_column(df: pd.DataFrame) -> pd.DataFrame:
        return file_processor_utils.add_reporter_intensity_corrected_1_column(df)

    @staticmethod
    def filter_out_rows_with_no_reporter_ions(df: pd.DataFrame) -> pd.DataFrame:
        return file_processor_utils.filter_out_rows_with_no_reporter_ions(df)