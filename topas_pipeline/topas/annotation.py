import pandas as pd

from .. import utils


def read_topas_annotations(topas_annotation_file: str) -> pd.DataFrame:
    topas_annotation_df = pd.read_excel(topas_annotation_file)
    topas_annotation_df = utils.whitespace_remover(topas_annotation_df)

    topas_annotation_df["TOPAS_subscore_level"] = (
        topas_annotation_df["TOPAS_SUBSCORE"] + " - " + topas_annotation_df["LEVEL"]
    )
    # empty cell in WEIGHT column means weight = 1
    topas_annotation_df["WEIGHT"] = topas_annotation_df["WEIGHT"].fillna(1)

    topas_annotation_df = topas_annotation_df.rename(
        {
            "GROUP": "group",
            "LEVEL": "level",
            "TOPAS_SCORE": "TOPAS_score",
            "TOPAS_SUBSCORE": "TOPAS_subscore",
            "SCORING RULE": "Scoring rule",
            "WEIGHT": "weight",
            "GENE NAME": "Gene names",
            "MODIFIED SEQUENCE": "Modified sequence",
        },
        axis=1,
    )

    topas_sheet_sanity_check(topas_annotation_df)

    return topas_annotation_df


def topas_sheet_sanity_check(topas_annotation_df: pd.DataFrame) -> None:
    """
    Explain

    TODO: add check that there are no protein groups in the gene names column, e.g. Gene1;Gene2
    """
    scoring_rules_found = set(topas_annotation_df["Scoring rule"].str.lower().unique())
    valid_scoring_rules = {
        "highest z-score",
        "highest z-score (p-site)",
        "highest protein phosphorylation score (2nd level z-score, fh)",
        "highest kinase score (2nd level z-score, fh)",
        "summed z-score",
    }
    unknown_scoring_rules = list(scoring_rules_found - valid_scoring_rules)
    if len(unknown_scoring_rules) > 0:
        raise ValueError(f"Unknown scoring rules: {unknown_scoring_rules}")

    # validate that scoring rule is "highest z-score (p-site)" if modified sequence column is not empty
    scoring_rules_found = set(
        topas_annotation_df[~topas_annotation_df["Modified sequence"].isnull()][
            "Scoring rule"
        ]
        .str.lower()
        .unique()
    )
    valid_scoring_rules = {"highest z-score (p-site)"}
    unknown_scoring_rules = list(scoring_rules_found - valid_scoring_rules)
    if len(unknown_scoring_rules) > 0:
        raise ValueError(
            f"Invalid scoring rule for entry with modified sequence: {unknown_scoring_rules}"
        )
