import os
import sys
import time
import argparse
import traceback
import logging

from . import __version__, __copyright__, __git_commit_hash__
from .utils import init_file_logger, send_slack_message
from . import config
from . import simsi
from . import preprocess
from . import clinical_annotation
from . import report_creation
from . import metrics
from .topas import phosphorylation
from .topas import topas
from . import portal_updater

logger = logging.getLogger(__name__)


def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-c", "--config", required=True, help="Absolute path to configuration file."
    )
    args = parser.parse_args(argv)

    configs = config.load(args.config)

    # Create results folder and save configurations
    os.makedirs(configs.results_folder, exist_ok=True)
    init_file_logger(configs.results_folder, "Pipeline_log.txt")

    logger.info(f"TOPAS-pipeline version {__version__} {__git_commit_hash__}")
    logger.info(f"{__copyright__}")
    logger.info(
        f'Issued command: {os.path.basename(__file__)} {" ".join(map(str, argv))}'
    )
    logger.info("Pipeline started")

    with open(os.path.join(configs.results_folder, "configs.json"), "w") as jsonFile:
        jsonFile.write(configs.asjson())

    t0 = time.time()

    try:
        start_time = time.time()
        # 0) process MaxQuant results with SIMSI (~10 hours)
        simsi.run_simsi(
            results_folder=configs.results_folder,
            search_result_folder=configs.preprocessing.raw_data_location,
            sample_annotation_file=configs.sample_annotation,
            raw_file_folders=configs.raw_file_folders,
            simsi_config=configs.simsi,
            data_types=configs.data_types,
        )
        logger.info("--- %.1f seconds --- simsi" % (time.time() - start_time))

        start_time = time.time()
        # 1) preprocess data (~1.5 hours, mostly slow because of MaxLFQ)
        preprocess.preprocess_raw(
            results_folder=configs.results_folder,
            sample_annotation_file=configs.sample_annotation,
            metadata_annotation=configs.metadata_annotation,
            run_simsi=configs.simsi.run_simsi,
            simsi_folder=configs.simsi.simsi_folder,
            preprocessing_config=configs.preprocessing,
            data_types=configs.data_types,
        )
        logger.info("--- %.1f seconds --- preprocessing" % (time.time() - start_time))

        start_time = time.time()
        # 2) clinical processing (~3 minutes)
        clinical_annotation.add_clinical_annotations(
            results_folder=configs.results_folder,
            debug=configs.preprocessing.debug,
            clinic_proc_config=configs.clinic_proc,
            data_types=configs.data_types,
        )
        logger.info(
            "--- %.1f seconds --- clinical processing" % (time.time() - start_time)
        )

        start_time = time.time()
        # 3) compute rank, z-score, fold change and p-value (<1 minute)
        metrics.compute_metrics(
            results_folder=configs.results_folder,
            debug=configs.preprocessing.debug,
            data_types=configs.data_types,
        )
        logger.info("--- %.1f seconds --- metrics" % (time.time() - start_time))

        start_time = time.time()
        # 4) Run WP2 scoring (<1 minute)
        phosphorylation.psite_scoring(
            results_folder=configs.results_folder,
            extra_kinase_annot=configs.clinic_proc.extra_kinase_annot,
            data_types=configs.data_types,
        )
        logger.info("--- %.1f seconds --- wp2 scoring" % (time.time() - start_time))

        start_time = time.time()
        # 5) compute TOPAS scores (<1 minute)
        topas.compute_topas_scores(
            results_folder=configs.results_folder,
            topas_annotation_file=configs.clinic_proc.prot_baskets,
            metadata_file=configs.metadata_annotation,
        )
        logger.info("--- %.1f seconds --- TOPAS scoring" % (time.time() - start_time))

        start_time = time.time()
        # 6) report creation (~18 minutes)
        report_creation.create_report(
            results_folder=configs.results_folder,
            debug=configs.preprocessing.debug,
            report_config=configs.report,
            topas_annotation_file=configs.clinic_proc.prot_baskets,
            data_types=configs.data_types,
        )
        logger.info("--- %.1f seconds --- report creation" % (time.time() - start_time))

        message = "Pipeline finished"
    except Exception as e:
        message = f"{type(e).__name__}: {e}"
        logger.info(message)
        logger.info(traceback.format_exc())

    send_slack_message(message, configs.results_folder, configs.slack)

    t1 = time.time()
    total = t1 - t0
    logger.info(f"Pipeline finished in {total:.1f} seconds")
    logger.info(f"Results have been written to {configs.results_folder}")

    try:
        portal_updater.main(
            results_folder=configs.results_folder,
            sample_annotation=configs.sample_annotation,
            metadata_annotation=configs.metadata_annotation,
            config_portal=configs.portal,
        )
    except Exception as e:
        message = f"{type(e).__name__}: {e}"
        logger.info(message)
        logger.info(traceback.format_exc())
        send_slack_message(message, configs.results_folder, configs.slack)


if __name__ == "__main__":
    main(sys.argv[1:])
