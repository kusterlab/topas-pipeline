import os
import sys
import time
import argparse
import traceback
import logging

from . import __version__, __copyright__, __git_commit_hash__
from .utils import init_file_logger, send_slack_message
from . import config
from .search_qc import maxquant_qc
from . import simsi
from .preprocess import preprocess
from . import clinical_annotation
from . import metrics
from .topas import protein_phosphorylation
from .topas import ck_substrate_phosphorylation
from .topas import rtk_substrate_phosphorylation
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
        # simsi.run_simsi(
        #     results_folder=configs.results_folder,
        #     search_result_folder=configs.preprocessing.raw_data_location,
        #     sample_annotation_file=configs.sample_annotation,
        #     raw_file_folders=configs.raw_file_folders,
        #     simsi_config=configs.simsi,
        #     data_types=configs.data_types,
        # )
        logger.info("--- %.1f seconds --- simsi" % (time.time() - start_time))

        start_time = time.time()
        # 0) generate QC statistics and plots for MaxQuant results
        maxquant_qc.maxquant_qc(
            results_folder=configs.results_folder,
            data_types=configs.data_types,
            sample_annot=configs.sample_annotation,
        )
        logger.info("--- %.1f seconds --- maxquant qc" % (time.time() - start_time))

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
        # 2) Run protein phosphorylation scoring (<1 minute)
        protein_phosphorylation.protein_phospho_scoring(
            results_folder=configs.results_folder,
            sample_annotation_file=configs.sample_annotation,
            metadata_file=configs.metadata_annotation,
        )
        logger.info(
            "--- %.1f seconds --- protein phospho scoring" % (time.time() - start_time)
        )

        start_time = time.time()
        # 3) add protein/p-site clinical annotations (~3 minutes)
        clinical_annotation.add_clinical_annotations(
            results_folder=configs.results_folder,
            clinic_proc_config=configs.clinic_proc,
            data_types=configs.data_types,
        )
        logger.info(
            "--- %.1f seconds --- clinical processing" % (time.time() - start_time)
        )

        start_time = time.time()
        # 4) compute rank, z-score, fold change and p-value (<1 minute)
        metrics.compute_metrics(
            results_folder=configs.results_folder,
            debug=configs.preprocessing.debug,
            data_types=configs.data_types,
        )
        logger.info("--- %.1f seconds --- metrics" % (time.time() - start_time))

        start_time = time.time()
        # 5) Run cytoplasmic kinase substrate phosphorylation scoring
        ck_substrate_phosphorylation.calculate_cytoplasmic_kinase_scores(
            results_folder=configs.results_folder,
            metadata_file=configs.metadata_annotation,
            topas_kinase_substrate_file=configs.clinic_proc.topas_kinase_substrate_file,
            expression_corrected_input=True,
        )
        logger.info(
            "--- %.1f seconds --- ck substrate phosphorylation scoring"
            % (time.time() - start_time)
        )

        start_time = time.time()
        # 6) Run receptor tyrosine kinase substrate phosphorylation scoring
        rtk_substrate_phosphorylation.calculate_rtk_scores(
            results_folder=configs.results_folder,
            metadata_file=configs.metadata_annotation,
            extra_kinase_annot=configs.clinic_proc.extra_kinase_annot,
            fasta_file=configs.preprocessing.fasta_file,
        )
        logger.info(
            "--- %.1f seconds --- rtk substrate phosphorylation scoring"
            % (time.time() - start_time)
        )

        start_time = time.time()
        # 7) compute TOPAS scores (<1 minute)
        topas.compute_topas_scores(
            results_folder=configs.results_folder,
            topas_annotation_file=configs.clinic_proc.prot_baskets,
            metadata_file=configs.metadata_annotation,
        )
        logger.info("--- %.1f seconds --- TOPAS scoring" % (time.time() - start_time))

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
