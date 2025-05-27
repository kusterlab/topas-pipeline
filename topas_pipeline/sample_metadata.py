import os
import shutil
from pathlib import Path

import pandas as pd

from . import utils

"""
Patient metadata annotation file is an excel file which contains a mapping of sample names 
to their respective metadata. In contrast to the sample annotation file, sample replicates
are not included in this file.

Columns (subject to change semi-regularly):
- Sample name
- Paper_pseudo_identifier
- Program
- Entity
- code_oncotree
- breadcrumb_oncotree
- Histologic subtype
- Histologic subtype, specifications
- ICD03 - Morpho ID
- Timeline
- Batch_No
- TMT_Channel
- Tumor cell content
- TCC_Bioinf_MASTER
- DNA degradation
- RNA QC issue
- Tissue constraints
- comment
- QC
- Input Material
- Localisation
- Portal
- Empfehlung MTB
- Follow up 
- Sex
- year of birth
- Discussed retrospective
- Discussed prospective MTB
- tissue_topology
- tissue_topology_specification
- Paper Entity
- HD name
- RNALibrarySize
- RNAMappingRate
- RNADuplicationRate

"""


def find_metadata_path(metadata_file: str):
    if os.path.isfile(metadata_file):
        return metadata_file

    metadata_file_outdated = metadata_file.replace(
        "Retrospective_MTBs_Evaluation",
        "Retrospective_MTBs_Evaluation/Metadata_excel_outdated",
    )
    if os.path.isfile(metadata_file_outdated):
        return metadata_file_outdated

    raise ValueError(f"Could not find Metadata file: {metadata_file}")


def copy_metadata_file(metadata_file: str, results_folder: str):
    metadata_file = find_metadata_path(metadata_file)
    shutil.copyfile(metadata_file, Path(results_folder) / Path(metadata_file).name)


def load(metadata_file: str):
    metadata_file = find_metadata_path(metadata_file)

    try:
        metadata_df = pd.read_excel(metadata_file)
    except ValueError as e:
        if "Value must be either numerical or a string containing a wildcard" in str(e):
            raise ValueError(
                f"Cannot open metadata file, check if you have active column filters applied. {metadata_file}"
            )
    except PermissionError:
        raise PermissionError(
            f"Cannot open metadata file, check if you have it open in Excel. {metadata_file}"
        )

    metadata_df = utils.whitespace_remover(metadata_df)
    return metadata_df


def filter_metadata_column_as_list(
    metadata_df: pd.DataFrame, column: str, minimum_patients_per_entity: int
) -> list[any]:
    counts = metadata_df[column].value_counts()
    filtered_counts = counts.where(lambda x: x >= minimum_patients_per_entity)
    filtered_counts = filtered_counts.dropna()
    return filtered_counts.index.to_list()


def filter_by_metadata_column(
    metadata_df: pd.DataFrame, column: str, minimum_patients_per_entity: int
) -> pd.DataFrame:
    filtered_entries = filter_metadata_column_as_list(
        metadata_df, column, minimum_patients_per_entity
    )
    return metadata_df[metadata_df[column].isin(filtered_entries)]
