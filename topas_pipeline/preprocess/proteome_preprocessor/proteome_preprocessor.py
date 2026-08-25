from dataclasses import asdict
import logging
import pandas as pd
from pathlib import Path
import os
import time
from typing import List, Dict

from topas_pipeline.config import Preprocessing, Simsi
from topas_pipeline import identification_metadata as id_meta
from topas_pipeline.preprocess import phosphopeptides
from topas_pipeline.sample_annotation import (
    load_sample_annotation,
    filter_sample_annotation,
    get_channel_to_sample_id_dict,
)
from topas_pipeline.preprocess import picked_group
from topas_pipeline.preprocess import preprocess_tools as prep
from topas_pipeline.preprocess import sample_mapping
from topas_pipeline.preprocess.proteome_preprocessor.utils import RESULT_FILE_GETTER_REGISTRY
from topas_pipeline.preprocess.quant_data_loader.quant_data_loader_factory import (
    QuantDataLoaderFactory,
)
from topas_pipeline import utils
from topas_pipeline import sample_annotation, sample_metadata

logger = logging.getLogger("topas_pipeline" + "." + Path(__file__).stem)


class ProteomePreprocessor:
    def __init__(
        self,
        results_folder: str,
        sample_annotation_file: str,
        metadata_annotation_file: str,
        data_types: List[str],
        quant_strategy: str,
        quant_file_formats: Dict[str, str],
        preprocessing_config: Preprocessing,
        simsi_config: Simsi,
    ):
        """
        Initialize the ProteomePreprocessor with the provided configuration
        Args:
            results_folder (str): The folder to store preprocessing results
            sample_annotation_file (str): The path to the sample annotation file
            metadata_annotation_file (str): The path to the metadata annotation file
            data_types (List[str]): The list of data types to preprocess (e.g., ["fp", "pp"])
            quant_strategy (str): The quantification strategy to use (e.g., "LFQ", "TMT")
            quant_file_formats (Dict[str, str]): The mapping of data types to quant file formats
                (e.g., {"fp": "evidence", "pp": "diann"})
            preprocessing_config (Preprocessing): The preprocessing configuration
            simsi_config (Simsi): The SIMSI configuration
        """
        self.results_folder = results_folder
        self.sample_annotation_file = sample_annotation_file
        self.metadata_annotation_file = metadata_annotation_file
        self.data_types = data_types
        self.quant_strategy = quant_strategy
        self.quant_file_formats = quant_file_formats
        self.preprocessing_config = asdict(preprocessing_config)
        self.simsi_config = asdict(simsi_config)
        self.preprocessing_func = {"fp": self.preprocess_fp, "pp": self.preprocess_pp}

    def preprocess(self, overwrite: bool = True):
        """Preprocess proteome data based on the provided configuration"""
        sample_metadata.copy_metadata_file(self.metadata_annotation_file, self.results_folder)
        sample_annotation.copy_sample_annotation_file(
            self.sample_annotation_file, self.results_folder
        )
        prep.check_annot(
            self.results_folder,
            self.sample_annotation_file,
            self.metadata_annotation_file,
        )
        # Load and filter sample annotation dataframe
        self.sample_annotation_df = load_sample_annotation(self.sample_annotation_file)
        self.sample_annotation_df = filter_sample_annotation(
            self.sample_annotation_df, remove_qc_failed=True, remove_replicates=True
        )
        self.sample_annotation_df = self.sample_annotation_df.reset_index()

        # Check if FP and PP result file paths exist for each sample
        quit_flag = False
        for data_type in self.data_types:
            quant_file_format = self.quant_file_formats[data_type]
            results_files = self.get_results_file_list(
                data_type, quant_file_format
            )
            # Check if files exist for path in result_files
            missing_files = [f for f in results_files if not os.path.exists(f)]
            if missing_files:
                logger.error(
                    f"Missing {data_type} result files: {missing_files}. Please check the raw data folder and sample annotation."
                )
                quit_flag = True
        if quit_flag:
            raise FileNotFoundError("One or more result files are missing. Preprocessing aborted.")

        # For each data_type run its corresponding preprocess func
        for data_type in self.data_types:
            if os.path.exists(
                os.path.join(self.results_folder, f"preprocessed_{data_type}.csv")
            ):
                if not overwrite:
                    logger.info(
                        f"Preprocessing {data_type} skipped - found files already processed"
                    )
                    continue
                logger.info(f"Found existing results but overwrite flag was set.")

            logger.info(f"Preprocessing {data_type} starts")
            self.preprocess_proteome(
                data_type, self.quant_file_formats[data_type], self.sample_annotation_df
            )

    def preprocess_proteome(
        self, data_type: str, quant_file_format: str, sample_annotation_df: pd.DataFrame
    ):
        """Preprocess proteome data based on the data_type"""
        # Get result file list based on quant_file_format
        if (
            results_files := prep.read_evidence_file_list_from_cache(
                Path(self.results_folder), data_type
            )
        ) is None:
            # Can be 3 different types of result files: (evidence.txt, report.parquet, combined_ion.tsv)
            results_files = self.get_results_file_list(data_type, quant_file_format)
            prep.write_evidence_file_list_to_cache(
                Path(self.results_folder), data_type, results_files
            )

        # Check if preprocessed2 is available for datatype
        # If true: read csv as df else: load sample data
        preprocessed2_file = os.path.join(
            self.results_folder, f"preprocessed_{data_type}2.csv"
        )
        if os.path.exists(preprocessed2_file):
            logger.info(
                f"Reusing previously generated results for {data_type}: {preprocessed2_file}"
            )
            df = pd.read_csv(preprocessed2_file)
        else:
            quant_data_loader_class = QuantDataLoaderFactory.get_loader(
                self.quant_strategy
            )
            quant_data_loader = quant_data_loader_class(
                results_files,
                quant_file_format,
                sample_annotation_df,
                self.preprocessing_config["normalize_to_reference"],
                self.preprocessing_config["debug"],
                {
                    "results_folder": self.results_folder,
                    "data_type": data_type,
                    "quant_strategy": self.quant_strategy,
                    "simsi_folder": self.simsi_config["simsi_folder"],
                },
            )
            df = quant_data_loader.load_and_normalize()

            # Apply preprocess funtion based on datatype
            df = self.preprocessing_func[data_type](df)

            # Write df to results folder as preprocessed2
            df.to_csv(preprocessed2_file, index=False, float_format="%.6g")

        # One more step for FP and PP seperately (Check with Matthew if it is TMT specific)
        if data_type == "fp":
            filtered_sample_annotation_file = os.path.join(
                self.results_folder, "sample_annot_filtered.csv"
            )
            if "Batch" not in sample_annotation_df.columns:
                sample_annotation_df["Batch"] = sample_annotation_df["Batch FP"]
            channel_to_sample_id_dict = get_channel_to_sample_id_dict(
                sample_annotation_df,
                filtered_sample_annotation_file,
                remove_qc_failed=True,
                remove_replicates=False,
            )

            index_cols = utils.get_index_cols(data_type)
            df = sample_mapping.rename_columns_with_sample_ids(
                df, channel_to_sample_id_dict, index_cols=index_cols
            )
            df = df.set_index(index_cols)
        elif data_type == "pp":
            df = phosphopeptides.group_phosphopeptides_and_normalize(
                results_folder=self.results_folder,
                sample_annotation_file=self.sample_annotation_file,
                run_lfq=self.preprocessing_config["run_lfq"],
            )

        # Write df to results folder as preprocessed
        df.reset_index().to_csv(
            os.path.join(self.results_folder, f"preprocessed_{data_type}.csv"),
            index=False,
            float_format="%.6g",
        )

    def preprocess_fp(self, df: pd.DataFrame) -> pd.DataFrame:
        """Preprocess Full Proteome data"""
        logger.info("Preprocess fp function")

        # Apply picked protein group on gene level and filter at 1% FDR
        df = picked_group.picked_protein_grouping(
            df,
            self.results_folder,
            self.preprocessing_config["picked_fdr"],
            self.preprocessing_config["fasta_file"],
            self.preprocessing_config["fdr_num_threads"],
        )

        # Filter out contaminants and reverse sequences
        df = prep.filter_data(df, data_type="fp")

        # create columns to store metadata about the identifications, e.g. imputed, detected in batch, single peptide id
        df = id_meta.create_metadata_columns(df)

        # Mark number of peptide identifications per sample
        df = id_meta.mark_num_peptides(df)

        # Mark proteins with peptide identifications out of range for that protein
        df = id_meta.mark_peptide_id_out_of_range(df)

        # Mark proteins detected in the batch but not in the sample
        df = id_meta.mark_detected_in_batch(df)

        # log10 transform intensities and turn missing values into nans
        df = prep.log_transform_intensities(df)

        # Mark proteins with quant outside of dynamic range in batch (too low compared to max)
        df = id_meta.mark_quant_out_of_range(df)
        return df

    def preprocess_pp(self, df: pd.DataFrame) -> pd.DataFrame:
        """Preprocess Phospho proteome data"""
        logger.info("Preprocess_pp function")
        save_debug_df = prep.get_save_debug_df_function(
            self.preprocessing_config["debug"], self.results_folder, "pp"
        )

        # Re-map gene names based on uniprot identifiers in a fasta file. This is necessary
        # because MaxQuant uses their own uniprot->gene mapping file that cannot be changed.
        df = picked_group.remap_gene_names(df, self.preprocessing_config["fasta_file"])

        # create columns to store metadata about the identifications, e.g. imputed, detected in batch, single peptide id
        df = id_meta.create_metadata_columns(df)

        if self.preprocessing_config["imputation"]:
            # Impute missing values within batches
            df.to_csv(
                os.path.join(
                    self.results_folder, "preprocessed_pp_before_imputation.csv"
                ),
                index=False,
                float_format="%.4g",
            )
            df = prep.impute_data(df)
            save_debug_df(df, "_after_imputation")

        # Filter out contaminants, reverse sequences and non-phospho peptides
        df = prep.filter_data(df, data_type="pp")

        # Aggregate p-peptide intensities across fractions and charge states
        df = prep.sum_peptide_intensities(
            df, self.preprocessing_config["debug"], self.preprocessing_config["run_lfq"]
        )
        save_debug_df(df, "_after_aggregation")

        # log10 transform intensities and turn missing values into nans
        df = prep.log_transform_intensities(df)

        # convert to wide format, i.e. each column is a patient with its peptide abundances
        df = prep.convert_long_to_wide_format(df, has_metadata_cols=True)

        # Mark proteins outside of dynamic range in batch (too low compared to max)
        df = id_meta.mark_quant_out_of_range(df)

        # I think solution is to save columns of transfer as separate file and only take to use for reports
        # at least for now
        # test with not running simsi
        if len(df) > 0 and df.columns.str.startswith("Transferred spectra count").any():
            df.loc[:, df.columns.str.startswith("Transferred spectra count")].to_csv(
                os.path.join(self.results_folder, "Transfer_metadata.csv")
            )
            df = df.drop(
                df.loc[
                    :, df.columns.str.startswith("Transferred spectra count")
                ].columns,
                axis=1,
            )
        return df

    def get_results_file_list(
        self, data_type: str, quant_file_format: str
    ) -> List[str]:
        """Get the list of result files based on the quant_file_format and data_type"""
        raw_data_folder = self.preprocessing_config["raw_data_location"]
        func = RESULT_FILE_GETTER_REGISTRY.get(quant_file_format)
        if func is None:
            raise ValueError(f"Unsupported quant_file_format: {quant_file_format}")
        return func(
            data_type, raw_data_folder, self.sample_annotation_df
        )


if __name__ == "__main__":
    import argparse

    from topas_pipeline import config

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-c", "--config", required=True, help="Absolute path to configuration file."
    )
    parser.add_argument(
        "-o",
        "--overwrite",
        action="store_true",
        help="Ignore existing results and recompute outputs.",
    )
    args = parser.parse_args()

    configs = config.load(args.config)
    processor = ProteomePreprocessor(
        configs.results_folder,
        configs.sample_annotation,
        configs.metadata_annotation,
        configs.data_types,
        configs.quant_strategy,
        configs.quant_file_formats,
        configs.preprocessing,
        configs.simsi,
    )
    start_time = time.time()
    processor.preprocess(args.overwrite)
    logger.info("--- %.1f seconds --- preprocessing" % (time.time() - start_time))
