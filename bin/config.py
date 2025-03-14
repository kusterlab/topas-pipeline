import os
import json

from pathlib import Path

from .sample_annotation import load_sample_annotation

def load(config_file):
    with open(Path(__file__).parent.parent / "config_defaults.json", "r") as inp:
        config_defaults = json.load(inp)

    with open(config_file, "r") as inp:
        configs = json.load(inp)
    
    configs = merge_dicts(config_defaults, configs)

    return configs


def merge_dicts(defaults, custom):
    for key, value in custom.items():
        if isinstance(value, dict) and key in defaults and isinstance(defaults[key], dict):
            defaults[key] = merge_dicts(defaults[key], value)  # Recursive merge
        else:
            defaults[key] = value  # Overwrite default value
    return defaults


def validate_sample_annot(configs):
    # Check for sample + metadata annotation errors before starting pipeline, currently not used
    # _ = prep.check_annot(configs["sample_annotation"], configs["metadata_annotation"], prep.in_metadata)
    sample_annot_df = load_sample_annotation(configs["sample_annotation"])
    sample_annot_df.to_csv(
        os.path.join(configs["results_folder"], "sample_annotation.csv")
    )

    # check_config(configs)



# def check_config(configs):
#
#     # check data types
#     for data_type in configs["data_types"]:
#         if data_type.upper() not in ['FP', 'PP']:
#             raise ValueError(f'Data type {data_type} not accepted. Accepted data types are `fp` (full proteome) and `pp` (phospho).')
#
#     # check existence of file locations
#     if not os.path.exists(configs["sample_annotation"]):
#         raise ValueError(f'Sample annotation file {configs["sample_annotation"]} does not exist in this location.')
#     metadata_file = configs["metadata_annotation"]
#     if not os.path.isfile(metadata_file):
#             metadata_file = metadata_file.replace('Retrospective_MTBs_Evaluation',
#                                                   'Retrospective_MTBs_Evaluation/Metadata_excel_outdated')
# if not os.path.isfile(metadata_file):
#     raise ValueError(f'Metadata: {metadata_file} not found.')
