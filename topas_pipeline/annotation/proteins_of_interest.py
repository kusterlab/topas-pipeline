import pandas as pd

from .. import utils


def load_poi_annotation_df(poi_annotation_file: str) -> pd.DataFrame:
    """
    Load all proteins of interest annotations as a DataFrame.
    """
    topas_annotations_df = pd.read_excel(poi_annotation_file)
    return topas_annotations_df


def merge_with_poi_annotations_inplace(
    df: pd.DataFrame, poi_annotation_df: pd.DataFrame
):
    utils.merge_by_delimited_field(
        df,
        poi_annotation_df[
            ["Gene names", "POI_REPORT", "POI_EXPLORATORY", "POI_PRODICT"]
        ],
        field_name="Gene names",
        inplace=True,
    )
