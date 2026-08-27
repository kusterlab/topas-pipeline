import argparse
import logging
import os
from pathlib import Path
import re
from typing import Dict, List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from tqdm import tqdm

from topas_pipeline import config as config
from topas_pipeline.sample_annotation import load_sample_annotation
from topas_pipeline.preprocess.proteome_preprocessor.utils import (
    RESULT_FILE_GETTER_REGISTRY,
)
from topas_pipeline.preprocess.file_processor.result_file_loader_factory import (
    ResultFileLoaderFactory,
)

logger = logging.getLogger("PhosDx_QC")
logging.basicConfig(level=logging.INFO)


def load_plate_data(
    data_type: str,
    file_format: str,
    plate_sample_annotation_df: pd.DataFrame,
    raw_data_folder: str,
):
    """
    Load Sample data of a plate for a given data type and file format.
    Args:
        data_type (str): The type of data (e.g., 'fp','pp').
        file_format (str): The format of the quantification file (e.g., 'diann','ionquant').
        plate_sample_annotation_df (pd.DataFrame): Sample annotation for the plate.
        raw_data_folder (str): Path to the folder containing raw data files.
    Returns:
        pd.DataFrame: Concatenated DataFrame containing the relevant data for the plate.
    """
    getter_func = RESULT_FILE_GETTER_REGISTRY[file_format]
    file_paths = getter_func(
        data_type, raw_data_folder, plate_sample_annotation_df
    )
    file_loader = ResultFileLoaderFactory.get_loader("lfq_" + file_format)
    all_batches = file_loader(
        file_paths, sample_annotation_df=plate_sample_annotation_df
    ).load()
    all_batches = [
        x[(x["Reverse"] != "+") & (x["Potential contaminant"] != "+")]
        for x in all_batches
    ]
    usecols = [
        "Modified sequence",
        "Batch",
        "Reporter intensity corrected 1",
        "Leading proteins",
        "Charge",
        "RT",
    ]
    concat_df = pd.concat([x[usecols] for x in all_batches])
    concat_df["Batch"] = (
        concat_df["Batch"].str.split("Batch").apply(lambda x: x[1])
    )
    return concat_df


def plot_bar_peptide_count(df_long: pd.DataFrame, title: str):
    """
    Plot a bar chart showing the count of unique peptides across batches.
    Args:
        df_long (pd.DataFrame): Long-format DataFrame containing 'Modified sequence', 'Batch', and 'Reporter intensity corrected 1' columns.
        title (str): Title of the plot.
    Returns:
        matplotlib.figure.Figure: The generated figure object.
    """
    pivot_df = df_long.pivot_table(
        index="Modified sequence",
        columns="Batch",
        values="Reporter intensity corrected 1",
        aggfunc="sum",
    )
    fig, ax = plt.subplots(figsize=(16, 10))
    pivot_df.count().plot.bar(ax=ax)
    ax.axhline(
        y=pivot_df.count().median(),
        color="k",
        linestyle="--",
        label="median of medians",
    )
    ax.set_ylabel("Count")
    ax.set_title(title)
    ax.tick_params(axis="x", rotation=90)
    ax.legend(loc="upper right")
    fig.tight_layout()
    return fig


def plot_bar_protein_count(df_long: pd.DataFrame, title: str):
    """
    Plot a bar chart showing the count of proteins across batches.
    Args:
        df_long (pd.DataFrame): Long-format DataFrame containing 'Leading proteins', 'Batch', and 'Reporter intensity corrected 1' columns.
        title (str): Title of the plot.
    Returns:
        matplotlib.figure.Figure: The generated figure object.
    """
    pivot_df = df_long.pivot_table(
        index="Leading proteins",
        columns="Batch",
        values="Reporter intensity corrected 1",
        aggfunc="sum",
    )
    fig, ax = plt.subplots(figsize=(16, 10))
    pivot_df.count().plot.bar(ax=ax)
    ax.axhline(
        y=pivot_df.count().median(),
        color="k",
        linestyle="--",
        label="median of medians",
    )
    ax.set_ylabel("Count")
    ax.set_title(title)
    ax.legend(loc="upper right")
    fig.tight_layout()
    return fig


def plot_bar_peptide_sum_intensity(df_long: pd.DataFrame, title: str):
    """
    Plot a bar chart showing the sum of peptide intensities across batches.
    Args:
        df_long (pd.DataFrame): Long-format DataFrame containing 'Modified sequence', 'Batch', and 'Reporter intensity corrected 1' columns.
        title (str): Title of the plot.
    Returns:
        matplotlib.figure.Figure: The generated figure object.
    """
    pivot_df = df_long.pivot_table(
        index="Modified sequence",
        columns="Batch",
        values="Reporter intensity corrected 1",
        aggfunc="sum",
    )
    fig, ax = plt.subplots(figsize=(16, 10))
    np.log10(pivot_df.sum()).plot.bar(ax=ax)
    ax.axhline(
        y=np.log10(pivot_df.sum().median()),
        color="k",
        linestyle="--",
        label="median of medians",
    )
    ax.set_ylabel("Summed intensity (log10)")
    ax.set_title(title)
    ax.legend(loc="upper right")
    fig.tight_layout()
    return fig


def plot_box_peptide_log_intensity(df_long: pd.DataFrame, title: str):
    """
    Plot a box plot showing the distribution of peptide log intensities across batches.
    Args:
        df_long (pd.DataFrame): Long-format DataFrame containing 'Modified sequence', 'Batch', and 'Reporter intensity corrected 1' columns.
        title (str): Title of the plot.
    Returns:
        matplotlib.figure.Figure: The generated figure object.
    """
    pivot_df = df_long.pivot_table(
        index="Modified sequence",
        columns="Batch",
        values="Reporter intensity corrected 1",
        aggfunc="sum",
    )
    fig, ax = plt.subplots(figsize=(16, 10))
    sns.boxplot(data=np.log10(pivot_df), color="#1f77b4", ax=ax)
    ax.axhline(
        y=np.log10(pivot_df).median().median(),
        color="k",
        linestyle="--",
        label="median of medians",
    )
    ax.set_ylabel("Intensity (log10)")
    ax.set_title(title)
    ax.tick_params(axis="x", rotation=90)
    ax.legend(loc="upper right")
    fig.tight_layout()
    return fig


def plot_stacked_bar_missed_cleavages(df_long: pd.DataFrame, title: str):
    """
    Plot a stacked bar chart showing the distribution of missed cleavages across batches.
    Args:
        df_long (pd.DataFrame): Long-format DataFrame containing 'Modified sequence', 'Batch', and 'Reporter intensity corrected 1' columns.
        title (str): Title of the plot.
    Returns:
        matplotlib.figure.Figure: The generated figure object.
    """
    pivot_df = df_long.pivot_table(
        index="Modified sequence",
        columns="Batch",
        values="Reporter intensity corrected 1",
        aggfunc="sum",
    )
    pivot_df["Missed_cleavages"] = pivot_df.index.map(
        lambda x: len(re.findall(r"[KR](?!P)", x[1:-2]))
    )
    grouped_df = pivot_df.groupby("Missed_cleavages").count()
    grouped_df = grouped_df / grouped_df.sum() * 100
    fig, ax = plt.subplots(figsize=(16, 10))
    grouped_df.T.plot(kind="bar", stacked=True, ax=ax)
    ax.set_ylabel("Percentage")
    ax.set_title(title)
    ax.legend(loc="upper right")
    fig.tight_layout()
    return fig


def plot_stacked_bar_precursor_charge(df_long: pd.DataFrame, title: str):
    """
    Plot a stacked bar chart showing the distribution of precursor charges across batches.
    Args:
        df_long (pd.DataFrame): Long-format DataFrame containing 'Modified sequence', 'Charge', 'Batch', and 'Reporter intensity corrected 1' columns.
        title (str): Title of the plot.
    Returns:
        matplotlib.figure.Figure: The generated figure object.
    """
    pivot_df = (
        df_long.pivot_table(
            index=["Modified sequence", "Charge"],
            columns="Batch",
            values="Reporter intensity corrected 1",
            aggfunc="sum",
        )
        .reset_index()
        .groupby("Charge")
        .count()
        .drop(columns=["Modified sequence"])
    )
    pivot_df = pivot_df / pivot_df.sum() * 100
    fig, ax = plt.subplots(figsize=(16, 10))
    pivot_df.T.plot(kind="bar", stacked=True, ax=ax)
    ax.set_ylabel("Percentage")
    ax.set_title(title)
    ax.legend(loc="upper right")
    fig.tight_layout()
    return fig


def plot_bar_enrichment_efficiency(df_long: pd.DataFrame, title: str):
    """
    Plot a bar chart showing the enrichment efficiency across batches.
    Args:
        df_long (pd.DataFrame): Long-format DataFrame containing 'Modified sequence', 'Batch', and 'Reporter intensity corrected 1' columns.
        title (str): Title of the plot.
    Returns:
        matplotlib.figure.Figure: The generated figure object.
    """
    pivot_df = df_long.pivot_table(
        index="Modified sequence",
        columns="Batch",
        values="Reporter intensity corrected 1",
        aggfunc="sum",
    )
    enriched_peptides = pivot_df[pivot_df.index.str.contains("p")]
    total_peptides = pivot_df
    enrichment_efficiency = (
        enriched_peptides.sum() / total_peptides.sum()
    ) * 100
    fig, ax = plt.subplots(figsize=(16, 10))
    enrichment_efficiency.plot.bar(ax=ax)
    ax.axhline(
        y=enrichment_efficiency.median(),
        color="k",
        linestyle="--",
        label="median of medians",
    )
    ax.set_ylabel("Enrichment Efficiency (%)")
    ax.set_title(title)
    ax.tick_params(axis="x", rotation=90)
    ax.legend(loc="upper right")
    fig.tight_layout()
    return fig


def plot_stacked_bar_quartile_sum_intensity(df_long: pd.DataFrame, title: str):
    """
    Plot a stacked bar chart showing the distribution of summed intensities across quartiles for each batch.
    Args:
        df_long (pd.DataFrame): Long-format DataFrame containing 'Modified sequence', 'Batch', and 'Reporter intensity corrected 1' columns.
        title (str): Title of the plot.
    Returns:
        matplotlib.figure.Figure: The generated figure object.
    """
    quartile_bins = pd.qcut(df_long["RT"], 4, labels=False)
    df_long["Quartile"] = quartile_bins
    pivot_df = df_long.pivot_table(
        index="Quartile",
        columns="Batch",
        values="Reporter intensity corrected 1",
        aggfunc="sum",
    )
    pivot_df_grouped = pivot_df / pivot_df.sum() * 100
    fig, ax = plt.subplots(figsize=(16, 10))
    pivot_df_grouped.T.plot(kind="bar", stacked=True, ax=ax)
    ax.set_ylabel("Summed Intensity (%)")
    ax.set_title(title)
    ax.legend(loc="upper right")
    fig.tight_layout()
    return fig


QC_PLOT_REGISTRY = {
    "peptide_count": (
        plot_bar_peptide_count,
        {"title": "Peptide ID Count in {data_type} Plate {plate_no}"},
    ),
    "peptide_sum_intensity": (
        plot_bar_peptide_sum_intensity,
        {"title": "Peptide Summed Intensity in {data_type} Plate {plate_no}"},
    ),
    "peptide_intensity": (
        plot_box_peptide_log_intensity,
        {"title": "Peptide Intensity in {data_type} Plate {plate_no}"},
    ),
    "protein_count": (
        plot_bar_protein_count,
        {"title": "Protein ID Count in {data_type} Plate {plate_no}"},
    ),
    "peptide_missed_cleavages": (
        plot_stacked_bar_missed_cleavages,
        {
            "title": "Peptide Missed Cleavages distribution in {data_type} Plate {plate_no}"
        },
    ),
    "precursor_charge": (
        plot_stacked_bar_precursor_charge,
        {
            "title": "Precursor Charge distribution in {data_type} Plate {plate_no}"
        },
    ),
    "phospho_enrichment_efficiency": (
        plot_bar_enrichment_efficiency,
        {
            "title": "Phospho Enrichment Efficiency (Intensity based) in {data_type} Plate {plate_no}"
        },
    ),
    "rt_quartile_sum_intensity": (
        plot_stacked_bar_quartile_sum_intensity,
        {
            "title": "Summed Intensity by RT Quartile in {data_type} Plate {plate_no}"
        },
    ),
}


class QCEngine:
    """
    A class to generate QC plots for different data types and plates based on the provided configuration.
    """

    def __init__(self, registry: Dict[str, tuple]):
        """
        Initialize the QCEngine with a registry of QC plot functions and their corresponding arguments.
        Args:
            registry (Dict[str, tuple]): A dictionary mapping QC plot names to their corresponding functions and arguments.
        """
        self.registry = registry

    def run_plate_qc(
        self,
        data_type: str,
        plate_sample_annotation_df: pd.DataFrame,
        quant_file_format: str,
        raw_data_folder: str,
        results_folder: str,
        qc_plots_list: List[str],
    ):
        """
        Generate QC plots for a specific plate and data type.
        Args:
            data_type (str): The type of data (e.g., 'fp','pp').
            plate_no (str): The plate number to process.
            plate_sample_annotation_df (pd.DataFrame): Sample annotation for the plate.
            quant_file_format (str): The format of the quantification file (e.g., 'diann','ionquant').
            raw_data_folder (str): Path to the folder containing raw data files.
            results_folder (str): Path to the folder where QC plots will be saved.
            qc_plots_list (List[str]): List of QC plot names to generate for the plate.
        """
        plate_no = plate_sample_annotation_df[
            f"Plate_No_{data_type.upper()}"
        ].iloc[0]
        concat_df = load_plate_data(
            data_type,
            quant_file_format,
            plate_sample_annotation_df,
            raw_data_folder,
        )
        for qc_plot in tqdm(qc_plots_list):
            plot_func, plot_kwargs = self.registry[qc_plot]
            plot_kwargs = (
                plot_kwargs.copy()
            )  # Create a copy to avoid modifying the original
            plot_kwargs["title"] = plot_kwargs["title"].format(
                plate_no=plate_no, data_type=data_type.upper()
            )
            fig = plot_func(concat_df, **plot_kwargs)
            fig.savefig(
                results_folder / f"{qc_plot}_{data_type}_plate_{plate_no}.png"
            )
            plt.close(fig)  # Close the figure to free up memory

    def run(
        self,
        qc_config: Dict[str, List[str]],
        results_folder: str,
        datatypes: List[str],
        quant_file_formats: Dict[str, str],
        raw_data_folder: str,
        sample_annotation_path: str,
    ):
        """
        Generate QC plots for the specified data types and plates based on the provided configuration.
        Args:
            qc_config (Dict[str, List[str]]): QC configuration specifying which plots to generate for each data type.
            results_folder (str): Path to the folder where QC plots will be saved.
            datatypes (List[str]): List of data types to process.
            quant_file_formats (Dict[str, str]): Dictionary mapping data types to their corresponding quantification file formats.
            raw_data_folder (str): Path to the folder containing raw data files.
            sample_annotation_path (str): Path to the sample annotation file.
        """
        results_folder = Path(results_folder) / "qc_plots"
        os.makedirs(results_folder, exist_ok=True)
        # Load the sample annotation file
        sample_annotation_df = load_sample_annotation(sample_annotation_path)
        sample_annotation_df = sample_annotation_df[:]
        for data_type in datatypes:
            logger.info(f"Processing data type: {data_type}")
            unique_plates = (
                sample_annotation_df[f"Plate_No_{data_type.upper()}"]
                .dropna()
                .unique()
            )
            logger.info(
                f"Found {len(unique_plates)} unique plates for data type {data_type}."
            )
            logger.info(f"Unique plates {data_type}: {unique_plates}")

            for plate_no in unique_plates:
                logger.info(
                    f"Processing plate number: {plate_no} for data type: {data_type}"
                )
                plate_sample_annotation_df = sample_annotation_df[
                    sample_annotation_df[f"Plate_No_{data_type.upper()}"]
                    == plate_no
                ]
                self.run_plate_qc(
                    data_type=data_type,
                    plate_sample_annotation_df=plate_sample_annotation_df,
                    quant_file_format=quant_file_formats[data_type],
                    raw_data_folder=raw_data_folder,
                    results_folder=results_folder,
                    qc_plots_list=qc_config[data_type],
                )
        logger.info(
            "QC plots generation completed for all data types and plates"
        )
        logger.info(f"QC plots saved in {results_folder}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-c",
        "--config",
        dest="config",
        required=False,
        help="Absolute path to configuration file.",
    )

    args = parser.parse_args()

    configs = config.load(args.config)
    QC_CONFIG = {
        "fp": [
            "peptide_count",
            "peptide_sum_intensity",
            "peptide_intensity",
            "protein_count",
            "peptide_missed_cleavages",
            "precursor_charge",
            "rt_quartile_sum_intensity",
        ],
        "pp": [
            "peptide_count",
            "peptide_sum_intensity",
            "peptide_intensity",
            "protein_count",
            "peptide_missed_cleavages",
            "precursor_charge",
            "phospho_enrichment_efficiency",
            "rt_quartile_sum_intensity",
        ],
    }
    qc_engine = QCEngine(QC_PLOT_REGISTRY)
    qc_engine.run(
        qc_config=QC_CONFIG,
        results_folder=configs.results_folder,
        datatypes=configs.data_types,
        quant_file_formats=configs.quant_file_formats,
        raw_data_folder=configs.preprocessing.raw_data_location,
        sample_annotation_path=configs.sample_annotation,
    )
