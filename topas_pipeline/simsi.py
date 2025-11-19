import io
import re
import os
import sys
import argparse
from pathlib import Path
import subprocess
import logging
from typing import List, Union
import zipfile

import pandas as pd
import simsi_transfer.main as simsi
from job_pool import JobPool

from . import __version__, __copyright__, __git_commit_hash__
import topas_pipeline.preprocess.preprocess_tools as prep
from . import meta_input_file as mi
from . import config
from .utils import init_file_logger

# hacky way to get the package logger instead of just __main__ when running as python -m topas_pipeline.simsi ...
logger = logging.getLogger(__package__ + "." + __file__)


def main(argv):

    # Parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-c",
        "--config",
        dest="config",
        required=True,
        help="Absolute path to configuration file.",
    )
    args = parser.parse_args(argv)

    configs = config.load(args.config)

    if not configs.simsi.run_simsi:
        logger.info(f"run_simsi flag is set to False, skipping SIMSI")
        return

    logger.info(f"TOPAS-pipeline-SIMSI version {__version__} {__git_commit_hash__}")
    logger.info(f"{__copyright__}")
    logger.info(
        f'Issued command: {os.path.basename(__file__)} {" ".join(map(str, argv))}'
    )

    os.makedirs(configs.results_folder, exist_ok=True)

    with open(f"{configs.results_folder}/configs.json", "w") as jsonFile:
        jsonFile.write(configs.asjson())

    run_simsi(
        results_folder=configs.results_folder,
        raw_data_location=configs.preprocessing.raw_data_location,
        sample_annotation=configs.sample_annotation,
        raw_file_folders=configs.raw_file_folders,
        simsi_config=configs.simsi,
        data_types=configs.data_types,
    )


def run_simsi(*args, **kwargs) -> None:
    # propagate logs from simsi_transfer package to our current log handlers
    logging.getLogger("simsi_transfer").handlers = logging.getLogger(
        __package__
    ).handlers

    data_types = kwargs.pop("data_types")
    processingPool = JobPool(
        processes=2, timeout=130000, write_progress_to_logger=True
    )  # 108,000 seconds = 30 hours
    for data_type in data_types:
        kwargs_with_data_type = kwargs.copy()
        kwargs_with_data_type["data_type"] = data_type

        processingPool.applyAsync(run_simsi_data_type, args, kwargs_with_data_type)
    processingPool.checkPool(printProgressEvery=1)


def run_simsi_data_type(
    results_folder: str,
    search_result_folder: str,
    sample_annotation_file: str,
    raw_file_folders: List[str],
    simsi_config: config.Simsi,
    data_type: str,
):
    if not simsi_config.run_simsi:
        logger.info(f"run_simsi flag is set to False, skipping SIMSI for {data_type}")
        return

    init_file_logger(results_folder, f"SIMSI_log_{data_type}.txt")

    logger.info(f"SIMSI started")

    summary_files = prep.get_summary_files(
        search_result_folder, data_type, sample_annotation_file
    )

    copy_raw_files(
        raw_file_folders, summary_files, data_type, simsi_config.simsi_folder
    )
    meta_input_file = mi.get_meta_input_file_path(results_folder, data_type)

    # manually create meta input files
    create_meta_input_file(
        summary_files,
        data_type,
        simsi_config.simsi_folder,
        meta_input_file,
        simsi_config.correction_factor_mapping_file,
    )

    meta_input_df = pd.read_csv(meta_input_file, sep="\t")

    result_folder_name = Path(results_folder).name
    simsi_output_folder = get_simsi_output_folder(
        simsi_config.simsi_folder, data_type, result_folder_name
    )
    simsi_cache_folder = get_simsi_cache_folder(simsi_config.simsi_folder, data_type)
    if matched_summaries_folder := find_matching_summaries_folder(
        simsi_cache_folder, meta_input_file
    ):
        logger.info(
            f"Found a summaries folder that matches the current list of folders: {matched_summaries_folder}"
        )
        logger.info("Skipping SIMSI processing.")
        return

    run_simsi_single(
        meta_input_file,
        simsi_output_folder,
        simsi_cache_folder,
        simsi_config.tmt_ms_level,
        simsi_config.stringencies,
        simsi_config.tmt_requantify,
        simsi_config.maximum_pep,
        simsi_config.num_threads,
    )

    store_results_for_reuse(
        simsi_output_folder, simsi_cache_folder, meta_input_file, result_folder_name
    )


def get_simsi_output_folder(simsi_folder: str, data_type: str, result_folder_name: str):
    """
    The SIMSI output folder contains the MaRaCluster results and a "summaries" folder
    which contains the SIMSI results in MaxQuant output format.
    """
    return (
        Path(simsi_folder)
        / Path("simsi_output")
        / Path(data_type.upper())
        / Path(result_folder_name)
    )


def get_simsi_cache_folder(simsi_folder: str, data_type: str):
    """
    The SIMSI cache folder stores the mzML, extracted TMT intensities and MaRaCluster
    dat files such that they can be reused in multiple SIMSI runs.
    """
    return Path(simsi_folder) / Path("simsi_output") / Path(data_type.upper())


def find_matching_summaries_folder(simsi_cache_folder: Path, meta_input_file: Path):
    """
    Iterates over folders starting with 'summaries_' to see if any of them contains results for the current list of folders
    """

    summary_folders_and_archives = list(simsi_cache_folder.glob("summaries_*"))
    summary_folders_and_archives.sort(key=os.path.getmtime, reverse=True)
    for s in summary_folders_and_archives:
        if s.is_dir():
            meta_input_file_other = s / Path(meta_input_file.name)
            if not (meta_input_file.is_file() and meta_input_file_other.is_file()):
                continue
        elif s.name.endswith(".zip"):
            try:
                archive = zipfile.ZipFile(s, "r")
            except zipfile.BadZipFile:
                logger.warning(f"Encountered corrupt zip file: {s}")

            file_path_in_archive = s.with_suffix("").name + "/" + meta_input_file.name
            if not file_path_in_archive in archive.namelist():
                continue
            meta_input_file_other = io.BytesIO(archive.read(file_path_in_archive))
            s = s.with_suffix("")
        else:
            continue

        if mi.meta_input_files_equal(meta_input_file, meta_input_file_other):
            return s

    return None


def get_correction_factor_files(
    experiments: List[str],
    correction_file_mapping_file: Path,
) -> List[str]:
    # has two columns: experiment, correction_factor_file
    correction_file_mapping_df = pd.read_csv(
        str(correction_file_mapping_file), sep="\t", index_col="experiment"
    )
    correction_files = map_experiments_to_files(
        experiments, correction_file_mapping_df, "correction_factor_file"
    )

    if None in correction_files:
        missing_experiments = get_missing_experiments(experiments, correction_files)
        raise ValueError(
            f"No correction factor file specified for experiments: {missing_experiments}"
        )

    return correction_files


def get_missing_experiments(experiments: List[str], correction_files: List[str]):
    return [
        experiment
        for experiment, correction_file in zip(experiments, correction_files)
        if correction_file is None
    ]


def map_experiments_to_files(
    experiments: List[str], correction_file_mapping_df: pd.DataFrame, key: str
):
    correction_file_mapping_df = drop_duplicate_indices(correction_file_mapping_df)
    experiment_mapping = correction_file_mapping_df.to_dict("index")
    return [
        experiment_mapping.get(experiment, {key: None})[key]
        for experiment in experiments
    ]


def drop_duplicate_indices(df):
    return df[~df.index.duplicated(keep="last")]


def find_simsi_evidence_file(
    results_folder: str, simsi_folder: str, data_type: str
) -> Union[str, io.BytesIO]:
    meta_input_file = mi.get_meta_input_file_path(results_folder, data_type.upper())
    simsi_cache_folder = get_simsi_cache_folder(simsi_folder, data_type.upper())
    simsi_summaries_folder = find_matching_summaries_folder(
        simsi_cache_folder, meta_input_file
    )
    if not simsi_summaries_folder:
        raise ValueError(
            "Could not find a SIMSI summaries folder with a matching list of input folders"
        )

    simsi_evidence_file = f"{str(simsi_summaries_folder)}/p10/p10_evidence.txt"
    logger.info(f"Found matching SIMSI summaries folder: {simsi_evidence_file}")

    if not os.path.isfile(simsi_evidence_file):
        simsi_evidence_file = extract_simsi_evidence_file_from_archive(
            simsi_evidence_file
        )

    return simsi_evidence_file


def extract_simsi_evidence_file_from_archive(simsi_evidence_file: Path) -> io.BytesIO:
    zip_file = Path(str(Path(simsi_evidence_file).parent.parent) + ".zip")
    if not zip_file.is_file():
        raise ValueError("Could not find SIMSI summaries folder nor archive.")

    logger.info(f"Extracting SIMSI evidence file from archive {zip_file}")
    archive = zipfile.ZipFile(zip_file, "r")

    file_path_in_archive = zip_file.with_suffix("").name + "/p10/p10_evidence.txt"
    if not file_path_in_archive in archive.namelist():
        raise ValueError("Could not find SIMSI evidence file in archive.")
    return io.BytesIO(archive.read(file_path_in_archive))


def store_results_for_reuse(
    simsi_output_folder: Path,
    simsi_cache_folder: Path,
    meta_input_file: Path,
    results_folder_name: str,
):
    """
    Rename the maracluster_output and summaries directories so we don't overwrite results
    """
    simsi_summaries_folder = simsi_output_folder / Path("summaries")
    simsi_summaries_folder_new = simsi_cache_folder / Path(
        f"summaries_{results_folder_name}"
    )
    simsi_summaries_folder.rename(simsi_summaries_folder_new)

    simsi_maracluster_folder = simsi_output_folder / Path("maracluster_output")
    simsi_maracluster_folder_new = simsi_cache_folder / Path(
        f"maracluster_output_{results_folder_name}"
    )
    simsi_maracluster_folder.rename(simsi_maracluster_folder_new)

    simsi_log_file = simsi_output_folder / Path("SIMSI.log")
    simsi_log_file_new = simsi_summaries_folder_new / Path("SIMSI.log")
    simsi_log_file.rename(simsi_log_file_new)

    # keep a copy in the summaries folder so we can check in future runs if we can re-use the results
    (simsi_summaries_folder_new / meta_input_file.name).write_text(
        meta_input_file.read_text()
    )

    # rmdir() will only delete if the directory is empty
    simsi_output_folder.rmdir()


def run_simsi_single(
    meta_input_file: Path,
    simsi_output_folder: Path,
    simsi_cache_folder: Path,
    tmt_ms_level: str,
    stringencies: int,
    tmt_requantify: bool,
    maximum_pep: int,
    num_threads: int,
):
    simsi_configs = [
        "--meta_input_file",
        str(meta_input_file),
        "--output_folder",
        str(simsi_output_folder),
        "--cache_folder",
        str(simsi_cache_folder),
        "--tmt_ms_level",
        str(tmt_ms_level),
        "--skip_annotated_clusters",
        "--skip_msmsscans",
        "--skip_msms",
        "--stringencies",
        str(stringencies),
        "--maximum_pep",
        str(maximum_pep),
        "--num_threads",
        str(num_threads),
        "--num_threads_per_precursor_bin",
        str(2),
    ]
    if tmt_requantify:
        simsi_configs.append("--tmt_requantify")
    simsi.main(simsi_configs)


def create_meta_input_file(
    summary_files: List[str],
    data_type: str,
    simsi_folder: str,
    meta_input_file: Path,
    correction_file_mapping_file: str,
):
    mq_txt_folders = [Path(f).parent for f in summary_files]
    experiments = []
    for summary in summary_files:
        summary_df = pd.read_csv(summary, sep="\t")
        experiments.append(summary_df["Experiment"][0])

    simsi_raw_folders = [
        get_simsi_raw_file_folder(simsi_folder, data_type, experiment)
        for experiment in experiments
    ]

    tmt_correction_files = None
    if len(correction_file_mapping_file) > 0:
        tmt_correction_files = get_correction_factor_files(
            experiments, correction_file_mapping_file
        )

    mi.write_meta_input_file(
        meta_input_file, mq_txt_folders, simsi_raw_folders, tmt_correction_files
    )


def copy_raw_files(raw_file_folders, summary_files, data_type, simsi_folder):
    """
    Copies raw files (skips already copied files) for all batches listed in sample_annotion based on the summary.txt of the MaxQuant search.
    """
    simsi_cache_folder = get_simsi_cache_folder(simsi_folder, data_type)
    for summary_file in summary_files:
        logger.debug(f"Copying raw files for {summary_file}")
        summary_df = pd.read_csv(summary_file, sep="\t")

        raw_files = summary_df["Raw file"]
        raw_files = raw_files.loc[raw_files != "Total"]

        experiment = summary_df["Experiment"][0]

        batch_simsi_folder = get_simsi_raw_file_folder(
            simsi_folder, data_type, experiment
        )
        if not os.path.exists(batch_simsi_folder):
            os.makedirs(batch_simsi_folder)

        for f in raw_files:
            # if raw file already exists as mzML don't copy
            if os.path.exists(os.path.join(simsi_cache_folder, "mzML", f"{f}.mzML")):
                continue

            success = False
            for raw_file_folder in raw_file_folders:
                success = copy_with_subprocess(
                    os.path.join(raw_file_folder, f"{f}.raw"),
                    os.path.join(batch_simsi_folder, f"{f}.raw"),
                )
                if success:
                    break

            if not success:
                raise RuntimeError(f"Failed to find or copy the raw file {f}.")


def get_simsi_raw_file_folder(simsi_folder, data_type, experiment):
    return os.path.join(simsi_folder, "raw_files", data_type.upper(), experiment)


def copy_with_subprocess(source, dest):
    if not os.path.exists(source):
        logger.debug(f"Could not find source file {source}")
        return False

    if os.path.exists(dest) and os.stat(source).st_size == os.stat(dest).st_size:
        # logger.debug(f"Skipping dest file {dest}, already copied")
        return True

    logger.debug(" ".join(["cp", source, dest]))
    with subprocess.Popen(
        ["cp", source, dest], stdout=subprocess.PIPE, stderr=subprocess.PIPE
    ) as proc:
        return proc.wait() == 0


"""
python3 -m topas_pipeline.simsi -c config_patients.json
"""
if __name__ == "__main__":
    main(sys.argv[1:])
