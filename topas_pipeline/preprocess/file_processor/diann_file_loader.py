from concurrent.futures import ProcessPoolExecutor
from functools import partial
import logging
import pandas as pd
import re
from typing import List, Union

from topas_pipeline.io.reader import ReaderFactory
from topas_pipeline.preprocess.file_processor.constants import (
    DIANN_MOD_DICT,
    PSM_TSV_MOD_DICT,
)
from topas_pipeline.preprocess.file_processor.result_file_loader import (
    BaseResultFileLoader,
)
from topas_pipeline.preprocess.file_processor.utils import (
    convert_mods_to_mq_format,
    convert_mod_column_to_mq_format,
)

logger = logging.getLogger(__name__)


class DiannFileLoader(BaseResultFileLoader):
    def __init__(self, result_files: Union[str, List[str]], **kwargs):
        self.result_files = result_files
        self.report_df_index_cols = ["Modified sequence", "Charge"]
        self.options = kwargs

    def load(self, usecols: List[str] = None) -> List[pd.DataFrame]:
        """Load DIANN report files and their corresponding PSM files in parallel and return a list of merged DataFrames"""
        with ProcessPoolExecutor(
            max_workers=8
        ) as executor:  # TODO: make max_workers configurable
            all_batches = executor.map(
                partial(self.load_single_file, usecols=usecols),
                self.result_files,
            )
        return list(all_batches)

    def load_single_file(
        self, diann_report_file_path: str, usecols: List[str] = None
    ) -> pd.DataFrame:
        """Load a single DIANN report file and its corresponding PSM file, merge them, and return the resulting DataFrame"""
        try:
            # Load the DIANN report file and convert it to evidence format
            report_df = self.load_diann_report_file(
                diann_report_file_path, usecols=usecols
            )

            # Get path to the corresponding PSM file based on the DIANN report file path
            psm_file_path = self.extract_psm_file_path_from_diann_report_path(
                diann_report_file_path
            )

            # Load the PSM file and convert it to evidence format
            psm_by_precursor_df = self.load_psm_file(psm_file_path)

            # Merge the PSM information into the evidence dataframe
            merged_df = pd.merge(
                report_df, psm_by_precursor_df, on=self.report_df_index_cols
            )

            # Additional processing for evidence dataframe
            merged_df = self.modify_protein_and_peptide_info(merged_df)

            # Update batch information
            merged_df = self.update_batch_info(
                merged_df, self.options["sample_annotation_df"].iloc[0]["Cohort"]
            )

            # Add reporter intensity corrected 1 column
            merged_df = self.add_reporter_intensity_corrected_1_column(merged_df)

            return merged_df
        except Exception as e:
            logger.error(f"Error processing file {diann_report_file_path}: {e}")
            return None

    def load_diann_report_file(
        self, diann_report_file_path: str, usecols: List[str] = None
    ) -> pd.DataFrame:
        """Load a DIANN report.tsv file and convert it to evidence format"""
        reader = ReaderFactory.get_reader(diann_report_file_path)
        report_df = reader.read(usecols=usecols)
        report_df = self.convert_diann_report_df_to_evidence_format(
            report_df, self.report_df_index_cols
        )
        report_df["Experiment"] = self.extract_experiment_name_from_diann_report_path(
            diann_report_file_path
        )
        return report_df

    @staticmethod
    def convert_diann_report_df_to_evidence_format(
        report_df: pd.DataFrame, index_cols: List[str]
    ) -> pd.DataFrame:
        """Convert a DIANN report DataFrame to evidence format"""
        report_df["Modified sequence"] = convert_mods_to_mq_format(
            report_df["Modified.Sequence"], DIANN_MOD_DICT
        )
        report_df["Charge"] = report_df["Precursor.Charge"]
        report_df["Raw file"] = report_df["Run"]
        report_df["id"] = report_df.index
        report_df["Intensity"] = report_df["Ms1.Normalised"]
        report_df = report_df[
            report_df["Intensity"] > 0.0
        ]  # up to 10% have 0 intensity, probably identified by MS/MS but no MS1 feature found

        # report_df["Experiment"] = sample_name # Add it after loading the report file
        report_df["Type"] = "MULTI-MSMS"
        report_df["Fraction"] = 1
        agg_funcs = {
            "PEP": "max", # TODO: check if this is the right aggregation function for PEP
            "Raw file": "first",
            "id": "first",
            "Intensity": "sum",
            # "Experiment": "first",
            "Type": "first",
            "Fraction": "first",
        }
        report_df = (
            report_df[list(agg_funcs.keys()) + index_cols]
            .groupby(index_cols, sort=False)
            .agg(agg_funcs)
        ).reset_index()
        return report_df

    @staticmethod
    def extract_experiment_name_from_diann_report_path(
        diann_report_file_path: str,
    ) -> str:
        """Extract the experiment name from the DIANN report file path"""
        # Extract the experiment name from the DIANN report file path
        # experiment_name = re.search(r"(?<=_P)\d+", diann_report_file_path).group()
        # Experiment name is in the format P<batch_name>_<other_info>_R<replicate_number>
        # We want to extract P<batch_name> and R<replicate_number> if present, and join them with a hyphen
        experiment_name = "-".join(
            item
            for item in re.search(
                r"_P(\d+)(?:_[A-Za-z0-9]+)(?:_R(\d))?", diann_report_file_path
            ).groups()
            if item is not None
        )
        return experiment_name

    def extract_psm_file_path_from_diann_report_path(
        self, diann_report_file_path: str
    ) -> str:
        """Extract the corresponding PSM file path from the DIANN report file path"""
        # Extract the PSM file path from the DIANN report file path
        psm_file_path = diann_report_file_path.split("dia-quant-output")[0] + "psm.tsv"
        return psm_file_path

    def load_psm_file(self, psm_file_path: str) -> pd.DataFrame:
        """Load a PSM file and convert it to evidence format"""
        reader = ReaderFactory.get_reader(psm_file_path)
        psm_df = reader.read()
        psm_by_precursor_df = self.convert_psm_df_to_evidence_format(
            psm_df, self.report_df_index_cols
        )
        return psm_by_precursor_df

    @staticmethod
    def convert_psm_df_to_evidence_format(
        psm_df: pd.DataFrame, index_cols: List[str]
    ) -> pd.DataFrame:
        """Convert a PSM DataFrame to evidence format"""
        psm_df["Modified Peptide"] = psm_df["Modified Peptide"].fillna(
            psm_df["Peptide"]
        )
        psm_df["Modified Peptide"] = convert_mods_to_mq_format(
            psm_df["Modified Peptide"], PSM_TSV_MOD_DICT
        )
        psm_df["Reverse"] = psm_df["Is Decoy"].apply(lambda x: "+" if x else "")
        psm_df["Potential contaminant"] = psm_df["Is Contaminant"].apply(
            lambda x: "+" if x else ""
        )
        # psm_df["PEP"] = 1.0 - psm_df["Probability"]
        psm_df["Score"] = psm_df["Hyperscore"]
        psm_df["MS/MS scan number"] = (
            psm_df["Spectrum"].str.split(".").str[-2].astype(int)
        )
        psm_df["Modifications"] = (
            psm_df["Assigned Modifications"]
            .fillna("")
            .apply(convert_mod_column_to_mq_format)
        )
        psm_df = psm_df.rename(columns={"Modified Peptide": "Modified sequence"})

        psm_df["Leading proteins"] = psm_df["Protein"]
        psm_df["Proteins"] = (
            psm_df["Leading proteins"]
            + ";"
            + psm_df["Mapped Proteins"].fillna("").str.split(", ").str.join(";")
        ).str.strip(";")
        psm_df["Gene names"] = (
            psm_df["Gene"].fillna("")
            + ";"
            + psm_df["Mapped Genes"].fillna("").str.split(", ").str.join(";")
        ).str.strip(";")

        agg_funcs = {
            "Reverse": "first",
            "Potential contaminant": "first",
            # "PEP": "max",
            "Score": "max",
            "MS/MS scan number": "first",
            "Modifications": "first",
            "Leading proteins": "first",
            "Proteins": "first",
            "Gene names": "first",
        }
        psm_by_precursor_df = (
            psm_df[list(agg_funcs.keys()) + index_cols]
            .groupby(index_cols, sort=False)
            .agg(agg_funcs)
        ).reset_index()

        return psm_by_precursor_df

    @staticmethod
    def update_batch_info(df: pd.DataFrame, cohort_name: str) -> pd.DataFrame:
        """Update batch information in the DataFrame based on the file path"""
        df["Batch"] = df["Experiment"].apply(lambda x: f"{cohort_name}_Batch{x}")
        return df
