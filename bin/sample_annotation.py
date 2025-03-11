import logging
from pathlib import Path
from typing import Dict, List, Optional, Union

import pandas as pd

from bin import utils


logger = logging.getLogger(__name__)


"""
Sample annotation file is a comma separated file contains a mapping of sample names 
to their respective batch, cohort and TMT channel.

Mandatory columns:
- Sample name
- Batch Name
- TMT Channel
- Cohort

Optional columns:
- QC
- Replicate
- Material issue
- Program
- Entity
- Histologic subtype
"""


def load_sample_annotation(sample_annotation_file: Union[str, Path]) -> pd.DataFrame:
    try:
        sample_annotation_df = pd.read_csv(sample_annotation_file, dtype={"Batch Name": str})
        sample_annotation_df["Sample name"] = sample_annotation_df[
        "Sample name"
    ].str.strip()
        sample_annotation_df = sample_annotation_df.set_index("Sample name")
    except PermissionError:
        raise PermissionError(
            f"Cannot open sample annotation file, check if you have it open in Excel. {sample_annotation_file}"
        )
    return sample_annotation_df


def filter_sample_annotation(
    sample_annotation_df: pd.DataFrame, remove_qc_failed: bool, remove_replicates: bool = False, remove_reference: bool = False
) -> pd.DataFrame:
    sample_annotation_df = sample_annotation_df.reset_index()
    sample_annotation_df["Sample name"] = sample_annotation_df[
        "Sample name"
    ].str.strip()
    sample_annotation_df = sample_annotation_df.set_index("Sample name")
    if "QC" in sample_annotation_df.columns and remove_qc_failed:
        sample_annotation_df = sample_annotation_df[
            sample_annotation_df["QC"].isin(["passed", "shaky"])
        ]
    if "Material issue" in sample_annotation_df.columns:
        sample_annotation_df = sample_annotation_df[
            sample_annotation_df["Material issue"] != "+"
        ]
    if "Replicate" in sample_annotation_df.columns and remove_replicates:
        sample_annotation_df = sample_annotation_df[
            sample_annotation_df["Replicate"] != "replicate"
        ]
    if "is_reference" in sample_annotation_df.columns and remove_reference:
        sample_annotation_df = sample_annotation_df[
            sample_annotation_df["is_reference"] != True
        ]
    return sample_annotation_df


def get_unique_batches(sample_annotation_df: pd.DataFrame) -> List[str]:
    return (
        sample_annotation_df[["Batch Name", "Cohort"]].drop_duplicates().values.tolist()
    )


def filter_samples_by_metadata(
    sample_annotation_df: pd.DataFrame,
    program: str,
    entity: str,
    histologic_subtype: str,
) -> pd.DataFrame:
    # check that only one of program, entity and histologic subtype is chosen
    counter = 0
    subset_types = [program, entity, histologic_subtype]
    for i, subset_type in enumerate(subset_types):
        if subset_type.upper() not in ["ALL", " ", ""]:
            counter += 1
            if counter > 1:
                ValueError(
                    "Only one of program, entity, histologic subtype can be used for subsetting before scoring."
                )
            chosen_subtype_count = i
            chosen_subtype = subset_type
    if counter == 0:
        logger.info("Processing data for all programs, entities and subtypes")
    else:
        subset_columns = ["Program", "Entity", "Histologic subtype"]

        # In case more than one program, entity or histologic subtype is chosen
        if "," in chosen_subtype or ";" in chosen_subtype:
            chosen_subtypes = utils.split_str_list(chosen_subtype)
            if chosen_subtype_count == 0:
                chosen_subtypes = [subtype.upper() for subtype in chosen_subtypes]
            sample_annotation_df = sample_annotation_df[
                sample_annotation_df[subset_columns[chosen_subtype_count]].isin(
                    chosen_subtypes
                )
            ]
        else:
            sample_annotation_df = sample_annotation_df[
                sample_annotation_df[subset_columns[chosen_subtype_count]].str.upper()
                == subset_types[chosen_subtype_count].upper()
            ]
    return sample_annotation_df


def get_channel_to_sample_id_dict(
    sample_annotation_df: pd.DataFrame,
    filtered_sample_annotation_file: Optional[str] = None,
    remove_qc_failed: bool = True,
    remove_replicates: bool = False,
) -> Dict[str, str]:
    """Returns a dictionary mapping tmt channels in batches (e.g. "Reporter intensity corrected 3 Sarcoma_Batch4") to sample identifiers

    Args:
        sample_annotation_df (pd.DataFrame): sample annotation dataframe
        filtered_sample_annotation_file (Optional[str], optional): filename to write the filtered sample annotation to
        remove_qc_failed (bool, optional): remove QC failed channels. Defaults to True.
        remove_replicates (bool, optional): remove patient replicate channels. Defaults to False.

    Returns:
        Dict[str, str]: dictionary from channel name to sample name
    """
    # in the future we might want to take info from metadata file and if so the following might be combined
    # eg. (material issue is in the same column as qc)
    sample_annotation_df = filter_sample_annotation(
        sample_annotation_df, remove_qc_failed, remove_replicates
    )

    def add_leading_zero(x):
        if x < 10:
            return f'0{x}'
        else:
            return str(x)

    def generate_channel_name(x):
        # return f"Reporter intensity corrected {x['TMT Channel']} {x['Cohort']}_Batch{add_leading_zero(x['Batch Name'])}"
        return f"Reporter intensity corrected {x['TMT Channel']} {x['Cohort']}_Batch{x['Batch Name']}"

    sample_annotation_df["channel"] = sample_annotation_df[
        ["TMT Channel", "Cohort", "Batch Name"]
    ].apply(generate_channel_name, axis=1)

    if filtered_sample_annotation_file:
        sample_annotation_df.to_csv(filtered_sample_annotation_file)
    
    return dict(
        zip(
            sample_annotation_df["channel"].tolist(),
            sample_annotation_df.index.tolist(),
        )
    )



