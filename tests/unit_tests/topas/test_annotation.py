import numpy as np
import pandas as pd
import pytest

from topas_pipeline.topas import annotation


class TestReadTopasAnnotation:
    def test_reads_and_processes_valid_excel_file(self, mocker):
        mock_excel_data = pd.DataFrame(
            {
                "TOPAS_SCORE": ["topas1  ", "topas2"],
                "TOPAS_SUBSCORE": ["sub1", "sub2"],
                "SCORING RULE": [
                    "highest z-score",
                    "highest kinase score (2nd level z-score, fh)",
                ],
                "GROUP": ["(R)TK", "(R)TK"],
                "LEVEL": ["expression  ", "kinase activity"],
                "WEIGHT": [0.5, 0.8],
                "GENE NAME": ["gene1", "gene2  "],
                "MODIFIED SEQUENCE": [np.nan, np.nan],
            }
        )
        mocker.patch("pandas.read_excel", return_value=mock_excel_data)

        result = annotation.read_topas_annotations("dummy_path.xlsx")

        pd.testing.assert_frame_equal(
            result,
            pd.DataFrame(
                {
                    "TOPAS_score": ["topas1", "topas2"],
                    "TOPAS_subscore": ["sub1", "sub2"],
                    "Scoring rule": [
                        "highest z-score",
                        "highest kinase score (2nd level z-score, fh)",
                    ],
                    "group": ["(R)TK", "(R)TK"],
                    "level": ["expression", "kinase activity"],
                    "weight": [0.5, 0.8],
                    "Gene names": ["gene1", "gene2"],
                    "Modified sequence": [np.nan, np.nan],
                    "TOPAS_subscore_level": [
                        "sub1 - expression",
                        "sub2 - kinase activity",
                    ],
                }
            ),
        )


class TestTopasSheetSanityCheck:
    def test_valid_scoring_rules(self):
        data = {
            "Scoring rule": [
                "highest z-score",
                "highest z-score (p-site)",
                "highest protein phosphorylation score (2nd level z-score, fh)",
                "highest kinase score (2nd level z-score, fh)",
                "summed z-score",
            ],
            "Modified sequence": [None, None, None, None, None],
        }
        df = pd.DataFrame(data)

        try:
            annotation.topas_sheet_sanity_check(df)
        except ValueError:
            pytest.fail("basket_sheet_sanity_check raised ValueError unexpectedly!")

    def test_unknown_scoring_rules(self):
        data = {
            "Scoring rule": [
                "highest z-score",
                "highest z-score (p-site)",
                "highest protein phosphorylation score (2nd level z-score, fh)",
                "highest kinase score (2nd level z-score, fh)",
                "summed z-score",
                "my_fancy_new_scoring_method",
            ],
            "Modified sequence": [None, None, None, None, None, None],
        }
        df = pd.DataFrame(data)

        with pytest.raises(
            ValueError,
            match=r"Unknown scoring rules: \['my_fancy_new_scoring_method'\]",
        ):
            annotation.topas_sheet_sanity_check(df)

    def test_psite_only_highest_z_score_rule(self):
        data = {
            "Scoring rule": ["highest z-score", "highest z-score (p-site)"],
            "Modified sequence": ["seq1", "seq2"],
        }
        df = pd.DataFrame(data)

        with pytest.raises(
            ValueError,
            match="Invalid scoring rule for entry with modified sequence:",
        ):
            annotation.topas_sheet_sanity_check(df)
