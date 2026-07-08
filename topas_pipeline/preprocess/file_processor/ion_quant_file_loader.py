from concurrent.futures import ProcessPoolExecutor
from functools import partial
from glob import glob
import logging
import pandas as pd
import re
from typing import List, Union

from topas_pipeline.io.reader import ReaderFactory
from topas_pipeline.preprocess.file_processor.constants import (
    COMBINED_ION_TSV_MOD_DICT,
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


class IonQuantFileLoader(BaseResultFileLoader):
    def __init__(self, result_files: Union[str, List[str]], **kwargs):
        self.result_files = result_files
        self.ionquant_df_index_cols = ["Modified sequence", "Charge"]
        self.options = kwargs

    def load(self, usecols: List[str] = None) -> List[pd.DataFrame]:
        """Load IonQuant combined_ion files and their corresponding PSM files in parallel and return a list of merged DataFrames"""
        with ProcessPoolExecutor(
            max_workers=8
        ) as executor:  # TODO: make max_workers configurable
            all_batches = executor.map(
                partial(self.load_single_file, usecols=usecols),
                self.result_files,
            )
        return list(all_batches)

    def load_single_file(
        self, ionquant_file_path: str, usecols: List[str] = None
    ) -> pd.DataFrame:
        """Load a single IonQuant combined_ion file and its corresponding PSM files, merge them, and return the resulting DataFrame"""
        try:
            combined_ion_df = self.load_ionquant_combined_ion_file(
                ionquant_file_path, usecols=usecols
            )

            # Get path to the corresponding PSM files based on the IonQuant combined ion file path
            psm_file_paths = self.extract_psm_file_paths_from_ionquant_file_path(
                ionquant_file_path
            )

            # Load the PSM files and convert them to evidence format
            psm_by_precursor_df = self.load_psm_files(psm_file_paths)

            # Merge the PSM information into the evidence dataframe
            merged_df = pd.merge(
                combined_ion_df, psm_by_precursor_df, on=self.ionquant_df_index_cols
            )

            # Additional processing for evidence dataframe
            merged_df = self.modify_protein_and_peptide_info(merged_df)

            # Update batch information
            merged_df = self.update_batch_info(
                merged_df, self.options["sample_annotation_df"].iloc[0]["Cohort"]
            )

            return merged_df
        except Exception as e:
            logger.error(f"Error processing file {ionquant_file_path}: {e}")
            return None

    def load_ionquant_combined_ion_file(
        self, ionquant_file_path: str, usecols: List[str] = None
    ) -> pd.DataFrame:
        """Load a single IonQuant combined_ion file, convert it to evidence format, and return the resulting DataFrame"""
        reader = ReaderFactory.get_reader(ionquant_file_path)
        df = reader.read(usecols=usecols)
        df = self.convert_ionquant_df_to_evidence_format(df)
        df["Experiment"] = self.extract_experiment_name_from_ionquant_path(
            ionquant_file_path
        )
        return df

    @staticmethod
    def convert_ionquant_df_to_evidence_format(
        ionquant_df: pd.DataFrame,
    ) -> pd.DataFrame:
        """Convert an IonQuant combined_ion dataframe to evidence format"""
        intensity_cols = ionquant_df.columns[
            ionquant_df.columns.str.endswith(" Intensity")
        ]
        if len(intensity_cols) > 2:
            raise ValueError(
                f"Expected at most two columns ending with ' Intensity', found {len(intensity_cols)} columns instead."
            )

        ionquant_df["Intensity"] = ionquant_df.loc[:, intensity_cols].mean(axis=1)
        raw_files = [x.split()[0] for x in intensity_cols]
        spectral_count_columns = [f"{x} Spectral Count" for x in raw_files]
        # combined_ion_df["Reporter intensity count 1"] = combined_ion_df.loc[:, spectral_count_columns].sum(axis=1)

        column_mapping = {
            "Modified Sequence": "Modified sequence",
            "Protein": "Leading proteins",
            "Mapped Proteins": "Proteins",
            "Gene": "Leading gene",
            "Mapped Genes": "Gene names",
            "Charge": "Charge",
            "Intensity": "Intensity",
            # "Reporter intensity corrected 1": "Reporter intensity corrected 1",
            # "Reporter intensity count 1": "Reporter intensity count 1",
            # intensity_col: "Reporter intensity corrected 1",
            # f"{raw_file} Spectral Count": "Reporter intensity count 1",
        }
        ionquant_df = ionquant_df.loc[:, column_mapping.keys()].rename(
            columns=column_mapping
        )
        ionquant_df["Raw file"] = ";".join(raw_files)
        ionquant_df["id"] = ionquant_df.index
        # ionquant_df["Reporter intensity 1"] = ionquant_df["Reporter intensity corrected 1"]
        # ionquant_df["Intensity"] = ionquant_df["Reporter intensity corrected 1"]
        ionquant_df = ionquant_df[
            ionquant_df["Intensity"] > 0.0
        ]  # up to 10% have 0 intensity, probably identified by MS/MS but no MS1 feature found

        # ionquant_df["Experiment"] = sample_name
        ionquant_df["Proteins"] = (
            ionquant_df["Leading proteins"]
            + ";"
            + ionquant_df["Proteins"].fillna("").str.split(", ").str.join(";")
        ).str.strip(";")
        ionquant_df["Gene names"] = (
            ionquant_df["Leading gene"].fillna("")
            + ";"
            + ionquant_df["Gene names"].fillna("").str.split(", ").str.join(";")
        ).str.strip(";")
        ionquant_df["Type"] = "MULTI-MSMS"

        ionquant_df["Modified sequence"] = convert_mods_to_mq_format(
            ionquant_df["Modified sequence"], COMBINED_ION_TSV_MOD_DICT
        )

        return ionquant_df

    @staticmethod
    def extract_experiment_name_from_ionquant_path(ionquant_file_path: str) -> str:
        """Extract the experiment name from the IonQuant combined_ion file path"""
        experiment_name = "-".join(
            item
            for item in re.search(
                r"_P(\d+)(?:_[A-Za-z0-9]+)(?:_R(\d))?", ionquant_file_path
            ).groups()
            if item is not None
        )
        return experiment_name

    @staticmethod
    def update_batch_info(df: pd.DataFrame, cohort_name: str) -> pd.DataFrame:
        """Update the batch information in the DataFrame based on the cohort name"""
        df["Batch"] = df["Experiment"].apply(lambda x: f"{cohort_name}_Batch{x}")
        return df

    @staticmethod
    def extract_psm_file_paths_from_ionquant_file_path(
        ionquant_file_path: str,
    ) -> List[str]:
        """Extract the paths to the corresponding PSM files based on the IonQuant combined_ion.tsv file path"""
        root_folder = ionquant_file_path.split("combined_ion.tsv")[0]
        psm_file_paths = glob(f"{root_folder}*/psm.tsv")
        return psm_file_paths

    def load_psm_files(self, psm_file_paths: List[str]) -> pd.DataFrame:
        """Load PSM files and convert them to evidence format"""
        psm_dfs = []
        for psm_file_path in psm_file_paths:
            reader = ReaderFactory.get_reader(psm_file_path)
            psm_dfs.append(reader.read())
        psm_df = pd.concat(psm_dfs, axis=0)
        psm_by_precursor_df = self.convert_psm_df_to_evidence_format(
            psm_df, self.ionquant_df_index_cols
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
        psm_df["PEP"] = 1.0 - psm_df["Probability"]
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

        agg_funcs = {
            "Reverse": "first",
            "Potential contaminant": "first",
            "PEP": "max",
            "Score": "max",
            "MS/MS scan number": "first",
            "Modifications": "first",
        }
        psm_by_precursor_df = (
            psm_df[list(agg_funcs.keys()) + index_cols]
            .groupby(index_cols, sort=False)
            .agg(agg_funcs)
        )

        return psm_by_precursor_df
