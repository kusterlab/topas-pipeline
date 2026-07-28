import logging
import pandas as pd
from typing import List, Dict

from topas_pipeline.preprocess.file_processor.result_file_loader_factory import (
    ResultFileLoaderFactory,
)
from topas_pipeline.preprocess.quant_data_loader.tmt_quant_data_loader import (
    TMTQuantDataLoader,
)
from topas_pipeline.sample_annotation import load_sample_annotation
from topas_pipeline import simsi

logger = logging.getLogger(__name__)

class SimsiTMTQuantDataLoader(TMTQuantDataLoader):
    def __init__(
        self,
        results_file_list: List[str],
        quant_file_format: str,
        sample_annotation_df: pd.DataFrame,
        normalize_to_reference: bool,
        debug: bool,
        metadata: Dict,
    ):
        self.results_file_list = results_file_list
        self.quant_file_format = quant_file_format
        self.sample_annotation_df = sample_annotation_df
        self.normalize_to_reference = normalize_to_reference
        self.debug = debug
        self.metadata = metadata
        self.results_folder = self.metadata["results_folder"]
        self.data_type = self.metadata["data_type"]
        self.quant_strategy = self.metadata["quant_strategy"]
        self.simsi_folder = self.metadata["simsi_folder"]
        self.simsi_evidence_file = self.get_simsi_evidence_file(
            self.results_folder, self.simsi_folder, self.data_type
        )

        # Get Result File loader based on quant file format used (Maxquant, Ionquant, Diann)
        result_file_loader_class = ResultFileLoaderFactory.get_loader(
            "_".join([self.quant_strategy, self.quant_file_format])
        )
        self.result_file_loader = result_file_loader_class(
            self.simsi_evidence_file, sample_annotation_df=self.sample_annotation_df
        )
    
    @staticmethod
    def get_simsi_evidence_file(results_folder: str, simsi_folder: str, data_type: str):
        simsi_evidence_file = simsi.find_simsi_evidence_file(
            results_folder, simsi_folder, data_type
        )
        return simsi_evidence_file
