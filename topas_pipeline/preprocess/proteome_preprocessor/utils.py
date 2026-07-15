from pathlib import Path
import pandas as pd
from typing import List


def get_maxquant_result_files(
    data_type: str, raw_data_folder: str, sample_annotation_df: pd.DataFrame
) -> List[str]:
    """
    Get the list of MaxQuant result files based on the
    data_type and Batch Name from sample_annotation_df
    Args:
        data_type (str): The type of data (e.g., "fp" or "pp")
        raw_data_folder (str): The path to the raw data folder
        sample_annotation_df (pd.DataFrame): The sample annotation DataFrame containing Batch Name information
    Returns:
        List[str]: A list of paths to the MaxQuant result files (evidence.txt)
    Raises:
        ValueError: If the number of matched directories for a batch is not equal to 1
    """
    result_file_paths = []
    for _, row in sample_annotation_df.iterrows():
        batch_name = row[f"Batch Name"]
        data_type_results_folder = Path(raw_data_folder)
        matched_paths = [
            str(p)
            for p in data_type_results_folder.glob(
                f"*{batch_name}_{data_type.upper()}*"
            )
            if p.is_dir()
        ]
        if len(matched_paths) == 1:
            result_file_paths.append(f"{matched_paths[0]}/combined/txt/evidence.txt")
        else:
            raise ValueError(
                f"Expected one directory for batch {batch_name} in {data_type_results_folder}, but found {len(matched_paths)}"
            )
    return result_file_paths


def get_diann_result_files(
    data_type: str, raw_data_folder: str, sample_annotation_df: pd.DataFrame
) -> List[str]:
    """
    Get the list of DIANN result files (report.parquet) based on the
    data_type and Batch Name from sample_annotation_df
    Args:
        data_type (str): The type of data (e.g., "fp" or "pp")
        raw_data_folder (str): The path to the raw data folder
        sample_annotation_df (pd.DataFrame): The sample annotation DataFrame containing Batch Name information
    Returns:
        List[str]: A list of paths to the DIANN result files (report.parquet)
    Raises:
        ValueError: If the number of matched directories for a batch is not equal to 1
    """
    result_file_paths = []
    for _, row in sample_annotation_df.copy().iterrows():
        batch_name = row[f"Batch Name {data_type.upper()}"]
        experiment_name = row["Experiment"]
        if "-" in batch_name:
            batch_name, replicate_number = batch_name.split("-")
        else:
            replicate_number = None
        data_type_results_folder = (
            Path(raw_data_folder) / experiment_name / data_type.upper()
        )
        matched_paths = [
            str(p)
            for p in data_type_results_folder.glob(f"*P{batch_name}*")
            if p.is_dir()
        ]
        if replicate_number:
            matched_paths = [p for p in matched_paths if f"R{replicate_number}" in p]
        if len(matched_paths) == 1:
            result_file_paths.append(
                f"{matched_paths[0]}/dia-quant-output/report.parquet"
            )
        else:
            raise ValueError(
                f"Expected one directory for batch {batch_name} in {data_type_results_folder}, but found {len(matched_paths)}"
            )
    return result_file_paths


def get_ionquant_result_files(
    data_type: str, raw_data_folder: str, sample_annotation_df: pd.DataFrame
) -> List[str]:
    """
    Get the list of IonQuant result files (combined_ion.tsv) based on the
    data_type and Batch Name from sample_annotation_df
    Args:
        data_type (str): The type of data (e.g., "fp" or "pp")
        raw_data_folder (str): The path to the raw data folder
        sample_annotation_df (pd.DataFrame): The sample annotation DataFrame containing Batch Name information
    Returns:
        List[str]: A list of paths to the IonQuant result files (combined_ion.tsv)
    Raises:
        ValueError: If the number of matched directories for a batch is not equal to 1
    """
    result_file_paths = []
    for _, row in sample_annotation_df.copy().iterrows():
        batch_name = row[f"Batch Name {data_type.upper()}"]
        experiment_name = row["Experiment"]
        if "-" in row[f"Batch Name {data_type.upper()}"]:
            batch_name, replicate_number = batch_name.split("-")
        else:
            replicate_number = None
        data_type_results_folder = (
            Path(raw_data_folder) / experiment_name / data_type.upper()
        )
        matched_paths = [
            str(p)
            for p in data_type_results_folder.glob(f"*P{batch_name}*")
            if p.is_dir()
        ]
        if replicate_number:
            matched_paths = [p for p in matched_paths if f"R{replicate_number}" in p]
        if len(matched_paths) == 1:
            result_file_paths.append(f"{matched_paths[0]}/combined_ion.tsv")
        else:
            raise ValueError(
                f"Expected one directory for batch {batch_name} in {data_type_results_folder}, but found {len(matched_paths)}"
            )
    return result_file_paths
