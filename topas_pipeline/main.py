import os
import sys
import json
import time
import argparse
import traceback
import logging

from . import __version__, __copyright__, __git_commit_hash__
from .utils import init_file_logger, send_slack_message
from . import config
from . import simsi
from . import preprocess
from . import clinical_process
from . import report_creation
from . import metrics
from . import TOPAS_psite_scoring
from . import basket_scoring
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
    os.makedirs(configs["results_folder"], exist_ok=True)
    init_file_logger(configs["results_folder"], "Pipeline_log.txt")

    logger.info(f"TOPAS-pipeline version {__version__} {__git_commit_hash__}")
    logger.info(f"{__copyright__}")
    logger.info(
        f'Issued command: {os.path.basename(__file__)} {" ".join(map(str, argv))}'
    )
    logger.info("Pipeline started")

    jsonString = json.dumps(configs, indent=4)
    with open(os.path.join(configs["results_folder"], "configs.json"), "w") as jsonFile:
        jsonFile.write(jsonString)

    t0 = time.time()

    try:
        # 0) process MaxQuant results with SIMSI
        simsi.run_simsi(
            configs["results_folder"],
            configs["preprocessing"]["raw_data_location"],
            configs["sample_annotation"],
            configs["raw_file_folders"],
            **configs["simsi"],
            data_types=configs["data_types"],
            slack_webhook_url=configs["slack_webhook_url"],
        )

        # 1) preprocess data (~1.5 hours, mostly slow because of MaxLFQ)
        preprocess.preprocess_raw(
            configs["results_folder"],
            configs["sample_annotation"],
            configs["metadata_annotation"],
            configs["simsi"]["simsi_folder"],
            **configs["preprocessing"],
            data_types=configs["data_types"],
        )

        start_time = time.time()
        # 2) clinical processing (~3 minutes)
        clinical_process.clinical_process(
            configs["results_folder"],
            configs["preprocessing"]["debug"],
            **configs["clinic_proc"],
            data_types=configs["data_types"],
        )
        logger.info(
            "--- %s seconds --- clinical processing" % (time.time() - start_time)
        )

        start_time = time.time()
        # 3) compute rank, z-score, fold change and p-value (<1 minute)
        metrics.compute_metrics(
            configs["results_folder"],
            configs["preprocessing"]["debug"],
            data_types=configs["data_types"],
        )
        logger.info("--- %s seconds --- metrics" % (time.time() - start_time))

        start_time = time.time()
        # 4) Run WP2 scoring (<1 minute)
        TOPAS_psite_scoring.psite_scoring(
            configs["results_folder"],
            configs["clinic_proc"]["extra_kinase_annot"],
            data_types=configs["data_types"],
        )
        logger.info("--- %s seconds --- wp2 scoring" % (time.time() - start_time))

        start_time = time.time()
        # 5) compute basket scores (<1 minute)
        basket_scoring.compute_TOPAS_scores(
            configs["results_folder"],
            configs["preprocessing"]["debug"],
            data_types=configs["data_types"],
            baskets_file=configs["clinic_proc"]["prot_baskets"],
            metadata_file=configs["metadata_annotation"],
        )
        logger.info("--- %s seconds --- basket scoring" % (time.time() - start_time))

        start_time = time.time()
        # 6) report creation (~18 minutes)
        report_creation.create_report(
            configs["results_folder"],
            configs["preprocessing"]["debug"],
            **configs["report"],
            annot_file=configs["clinic_proc"]["prot_baskets"],
            data_types=configs["data_types"],
        )
        logger.info("--- %s seconds --- report creation" % (time.time() - start_time))

        message = "Pipeline finished"
    except Exception as e:
        message = f"{type(e).__name__}: {e}"
        logger.info(message)
        logger.info(traceback.format_exc())

    send_slack_message(message, configs["results_folder"], configs["slack_webhook_url"])

    t1 = time.time()
    total = t1 - t0
    logger.info(f"Pipeline finished in {total} seconds")
    logger.info(f'Results have been written to {configs["results_folder"]}')

    # Check if configs['portal']['update'] true and if so try to push results to portal
    try:
        portal_updater.main(configs)
    except Exception as e:
        message = f"{type(e).__name__}: {e}"
        logger.info(message)
        logger.info(traceback.format_exc())
        send_slack_message(
            message, configs["results_folder"], configs["slack_webhook_url"]
        )

    # After preprocess and clinical process: data stats, graphs etc
    # stats.analyse_results(**configs["results"])


if __name__ == "__main__":
    main(sys.argv[1:])
