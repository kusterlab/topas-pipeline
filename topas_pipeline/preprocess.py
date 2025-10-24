import os
import sys
import logging

import pandas as pd
from typing import Union
from pathlib import Path

from . import config
from . import preprocess_tools as prep
from . import picked_group as picked
from . import sample_annotation
from . import sample_metadata
from . import utils
from . import identification_metadata as id_meta
from .data_loaders.tmt_loader import TMTLoader
from .data_loaders.simsi_tmt_loader import SimsiTMTLoader
from .data_loaders.lfq_loader import LFQLoader

logger = logging.getLogger(__name__)


def preprocess_raw(**kwargs) -> None:
    # Unpack directories and do checks on annotation (sample annot, metadata)
    sample_annotation_file = kwargs.pop("sample_annotation_file")
    metadata_annotation = kwargs.pop("metadata_annotation")

    sample_metadata.copy_metadata_file(metadata_annotation, kwargs["results_folder"])
    sample_annotation.copy_sample_annotation_file(sample_annotation_file, kwargs["results_folder"])

    sample_annotation_df = prep.check_annot(kwargs["results_folder"],
        sample_annotation_file, metadata_annotation, prep.in_metadata
    )  # just for our pipeline

    data_types = kwargs.pop("data_types")
    for data_type in data_types:
        preprocess_raw_data_type(
            **kwargs,
            sample_annotation_df=sample_annotation_df,
            data_type=data_type,
        )


def preprocess_raw_data_type(
    results_folder: Union[str, Path],
    run_simsi: bool,
    simsi_folder: Union[str, Path],
    sample_annotation_df: pd.DataFrame,
    preprocessing_config: config.Preprocessing,
    data_type: str,
) -> None:
    """
    Calling function for preprocessing of both phospho and full proteome data
    :param results_folder: location in which to save preprocessed files
    :param run_simsi: boolean for whether to use simsi for identification transfer between batches
    :param simsi_folder: location in which simsi saves its processed files
    :param sample_annotation_df: df containing possible subfolder (of raw folder), metadata and sample information (tmt: batch and channel)
    :param raw_data_location: location of data to preprocess (MaxQuant search folders)
    :param picked_fdr: FDR value for picked group filtering
    :param fasta_file: used for protein grouping
    :param fdr_num_threads: number of threads to use for MaxLFQ tool
    :param program: subsetting of data for <program> during rest of pipeline when multiple is used for simsi and normalization
    :param entity: subsetting of data for <entity> during rest of pipeline when multiple is used for simsi and normalization
    :param histologic_subtype: subsetting of data for <histologic subtype> during rest of pipeline when multiple is used for simsi and normalization
    :param imputation: boolean for whether to impute or not (only for TMT phospho)
    :param run_lfq: boolean for when data is lfq else tmt is expected
    :param debug: boolean for saving more intermediate results files for debugging purpose
    :param data_type: full protome (fp) or phospho proteome (pp) for which data is analyzed (for lfq only fp)
    :param normalize_to_reference: boolean for whether or not to normalize intensities relative to the reference channels (e.g. if MS1 chromatography peaks are bad)
    """
    if not preprocessing_config.run_preprocessing:
        logger.info(
            f"run_preprocessing flag is set to False, skipping preprocessing for {data_type}"
        )
        return

    # check if file is there - if so skip this
    if os.path.exists(os.path.join(results_folder, f"preprocessed_{data_type}.csv")):
        logger.info(
            f"Preprocessing {data_type} skipped - found files already preprocessed"
        )
        return

    logger.info(f"Preprocessing {data_type} starts")

    preprocessed2_file = os.path.join(results_folder, f"preprocessed_{data_type}2.csv")
    if os.path.exists(preprocessed2_file):
        logger.info(
            f"Reusing previously generated results for {data_type}: {preprocessed2_file}"
        )
        df = pd.read_csv(preprocessed2_file)
    else:
        df = load_sample_data(
            results_folder,
            sample_annotation_df,
            run_simsi,
            simsi_folder,
            preprocessing_config.raw_data_location,
            preprocessing_config.run_lfq,
            preprocessing_config.debug,
            data_type,
            normalize_to_reference=preprocessing_config.normalize_to_reference,
        )

        preprocess_function = preprocess_fp
        if data_type == "pp":
            preprocess_function = preprocess_pp

        # ~5 minutes for pp, ~1 hour for fp, of which 1 hour is LFQ
        # returns dataframe in "wide format", i.e. patients as columns
        df = preprocess_function(
            df,
            results_folder,
            preprocessing_config.picked_fdr,
            preprocessing_config.fasta_file,
            preprocessing_config.fdr_num_threads,
            imputation=preprocessing_config.imputation,
            debug=preprocessing_config.debug,
            run_lfq=preprocessing_config.run_lfq,
        )
        df.to_csv(preprocessed2_file, index=False, float_format="%.6g")

    filtered_sample_annotation_file = os.path.join(
        results_folder, "sample_annot_filtered.csv"
    )
    channel_to_sample_id_dict = sample_annotation.get_channel_to_sample_id_dict(
        sample_annotation_df,
        filtered_sample_annotation_file,
        remove_qc_failed=True,
        remove_replicates=False,
    )

    index_cols = utils.get_index_cols(data_type)
    df = prep.rename_columns_with_sample_ids(
        df, channel_to_sample_id_dict, index_cols=index_cols
    )

    df = df.set_index("Gene names")

    df.reset_index().to_csv(
        os.path.join(results_folder, f"preprocessed_{data_type}.csv"),
        index=False,
        float_format="%.6g",
    )


def load_sample_data(
    results_folder: Union[str, Path],
    sample_annotation_df: pd.DataFrame,
    run_simsi: bool,
    simsi_folder: Union[str, Path],
    raw_data_location: Union[str, Path],
    run_lfq: bool,
    debug: bool,
    data_type: str,
    normalize_to_reference: bool,
) -> pd.DataFrame:
    evidence_files = prep.get_evidence_files(
        sample_annotation_df, raw_data_location, data_type
    )
    data_loader = TMTLoader(evidence_files)
    if run_simsi:
        data_loader = SimsiTMTLoader(
            evidence_files, results_folder, simsi_folder, data_type
        )
    elif run_lfq:
        data_loader = LFQLoader(evidence_files)

    df = prep.load_and_normalize(
        data_loader,
        results_folder,
        sample_annotation_df,
        data_type=data_type,
        normalize_to_reference=normalize_to_reference,
        debug=debug,
    )
    return df


def preprocess_pp(
    df: pd.DataFrame,
    results_folder: Union[str, Path],
    picked_fdr: float,
    fasta_file: str,
    fdr_num_threads: int,
    imputation: bool,
    debug: bool,
    run_lfq: bool,
) -> pd.DataFrame:
    logger.info("Preprocess_pp function")

    # Re-map gene names based on uniprot identifiers in a fasta file. This is necessary
    # because MaxQuant uses their own uniprot->gene mapping file that cannot be changed.
    df = picked.remap_gene_names(df, fasta_file)

    # create columns to store metadata about the identifications, e.g. imputed, detected in batch, single peptide id
    df = id_meta.create_metadata_columns(df)

    if imputation:
        # Impute missing values within batches
        df.to_csv(
            os.path.join(results_folder, "preprocessed_pp_before_imputation.csv"),
            index=False,
            float_format="%.4g",
        )
        df = prep.impute_data(df)
        if debug:
            df.to_csv(
                os.path.join(
                    results_folder, "debug_preprocessed_pp_after_imputation.csv"
                ),
                index=False,
                float_format="%.4g",
            )

    # Filter out contaminants, reverse sequences and non-phospho peptides
    df = prep.filter_data(df, data_type="pp")

    # Aggregate p-peptide intensities across fractions and charge states
    df = prep.sum_peptide_intensities(df, debug, run_lfq)
    if debug:
        df.to_csv(
            os.path.join(results_folder, "debug_preprocessed_pp_after_aggregation.csv"),
            index=False,
            float_format="%.4g",
        )

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
            os.path.join(results_folder, "Transfer_metadata.csv")
        )
        df = df.drop(
            df.loc[:, df.columns.str.startswith("Transferred spectra count")].columns,
            axis=1,
        )
    return df


def preprocess_fp(
    df: pd.DataFrame,
    results_folder: Union[str, Path],
    picked_fdr: float,
    fasta_file: str,
    fdr_num_threads: int,
    imputation: bool,
    debug: bool,
    run_lfq: bool,
) -> pd.DataFrame:
    """
    Function for preprocessing full proteome MQ evidence files or simsi results
    :param results_folder: location to save results from picked group fdr + log file
    :param picked_fdr: FDR value for filtering after picked group FDR
    :return: Dataframe with combined preprocessed data from all batches
    """
    logger.info("Preprocess fp function")

    # Apply picked protein group on gene level and filter at 1% FDR
    df = picked.picked_protein_grouping(
        df, results_folder, picked_fdr, fasta_file, fdr_num_threads
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

    # Mark proteins with quant outside of dynamic range in batch (too low compared to max)
    df = id_meta.mark_quant_out_of_range(df)

    # log10 transform intensities and turn missing values into nans
    df = prep.log_transform_intensities(df)
    return df


if __name__ == "__main__":
    import argparse

    from . import config

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-c", "--config", required=True, help="Absolute path to configuration file."
    )
    args = parser.parse_args(sys.argv[1:])

    configs = config.load(args.config)

    preprocess_raw(
        results_folder=configs.results_folder,
        sample_annotation_file=configs.sample_annotation,
        metadata_annotation=configs.metadata_annotation,
        run_simsi=configs.simsi.run_simsi,
        simsi_folder=configs.simsi.simsi_folder,
        preprocessing_config=configs.preprocessing,
        data_types=configs.data_types,
    )
