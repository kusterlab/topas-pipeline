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
