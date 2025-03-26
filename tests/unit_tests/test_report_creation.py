import numpy as np
import pandas as pd
import pytest

from topas_pipeline import report_creation


class TestMergeTopasScoreWithTopasSubscore:
    # merge TOPAS score and subscore names correctly when both are non-empty
    def test_merge_topas_score_and_subscore_names(self):
        import pandas as pd

        data = {"TOPAS_score": "fruit;vegetable", "TOPAS_subscore": "apple;carrot"}
        row = pd.Series(data)

        result = report_creation.merge_topas_score_and_subscore_names(row)

        assert result == "fruit - apple;vegetable - carrot"


class TestGetUniqueTopasNames:
    def test_handles_list_input_correctly(self):
        topas_score_names = ["apple", "banana", "apple", "orange"]
        expected_output = "apple;banana;orange"
        assert (
            report_creation.get_unique_topas_names(topas_score_names) == expected_output
        )

    def test_handles_string_input_correctly(self):
        topas_score_names = "apple;banana;apple;orange"
        expected_output = "apple;banana;orange"
        assert (
            report_creation.get_unique_topas_names(topas_score_names) == expected_output
        )

    # Handles float input correctly
    def test_handles_float_input_correctly(self):
        topas_score = 3.5
        expected_output = 3.5
        assert report_creation.get_unique_topas_names(topas_score) == expected_output
