import os
import sys
import logging

import pandas as pd
from pathlib import Path
from typing import Union

from . import config
from . import clinical_tools
from . import utils
from .annotation import proteins_of_interest as poi

logger = logging.getLogger(__name__)

TOPAS_SCORE_COLUMNS = clinical_tools.TOPAS_SCORE_COLUMNS
TOPAS_SUBSCORE_COLUMNS = clinical_tools.TOPAS_SUBSCORE_COLUMNS


def add_clinical_annotations(*args, **kwargs) -> None:
    data_types = kwargs.pop("data_types")

    for data_type in data_types:
        add_clinical_annotations_data_type(*args, **kwargs, data_type=data_type)


def add_clinical_annotations_data_type(
    results_folder: Union[str, Path],
    clinic_proc_config: config.ClinicProc,
    data_type: str,
) -> None:
    """
    Opens preprocessed data, annotates phospho and TOPAS scores

    :param results_folder: path to results folder to which preprocessed data can be found and clinically processed data saved
    :param debug:
    :param clinic_proc_config: paths to file used for annotating to phospho and TOPAS scpres
    :param data_type: 'fp' for full proteome, 'pp' for phospho proteome
    """
    if os.path.exists(os.path.join(results_folder, f"annot_{data_type}.csv")):
        logger.info(
            f"Clinical processing {data_type} skipped - found files already preprocessed"
        )
        return

    if data_type.startswith("fp"):
        index_col = "Gene names"
        keep_default_na = True
    else:
        index_col = "Modified sequence group"
        keep_default_na = False
    preprocessed_df = pd.read_csv(
        os.path.join(results_folder, f"preprocessed_{data_type}.csv"),
        index_col=index_col,
        keep_default_na=keep_default_na,
    )

    if data_type == "pp":
        preprocessed_df["Modified sequence"] = preprocessed_df.index.str.split(";")
        preprocessed_df = preprocessed_df.explode("Modified sequence")

        # ~30 minutes for 3000 samples
        logger.info("Annotating phospho sites")
        preprocessed_df = clinical_tools.add_phospho_annotations(
            preprocessed_df,
            clinic_proc_config,
        )

    annot_levels = list(TOPAS_SCORE_COLUMNS.keys()) + list(
        TOPAS_SUBSCORE_COLUMNS.keys()
    )

    # ~20 minutes for 3000 samples for phospho
    for annot_type in annot_levels:
        preprocessed_df = clinical_tools.add_topas_annotations(
            preprocessed_df,
            clinic_proc_config.prot_baskets,
            data_type=data_type,
            annot_type=annot_type,
        )
    
    if len(clinic_proc_config.proteins_of_interest_file) > 0:
        poi_annotation_df = poi.load_poi_annotation_df(clinic_proc_config.proteins_of_interest_file)
        poi.merge_with_poi_annotations_inplace(preprocessed_df, poi_annotation_df)

    preprocessed_df.to_csv(
        os.path.join(results_folder, f"annot_{data_type}.csv"), float_format="%.6g"
    )


def read_annotated_expression_file(
    results_folder: Union[str, Path],
    data_type: str,
) -> pd.DataFrame:
    return pd.read_csv(
        os.path.join(results_folder, f"annot_{data_type}.csv"),
        index_col=utils.get_index_cols(data_type),
    )


"""
python3 -m topas_pipeline.clinical_annotation -c config_patients.json
"""
if __name__ == "__main__":
    import argparse

    from . import config

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-c", "--config", required=True, help="Absolute path to configuration file."
    )
    args = parser.parse_args(sys.argv[1:])

    configs = config.load(args.config)

    add_clinical_annotations(
        results_folder=configs.results_folder,
        clinic_proc_config=configs.clinic_proc,
        data_types=configs.data_types,
    )
