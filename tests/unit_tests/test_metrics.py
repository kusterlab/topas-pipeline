import os
from pathlib import Path

import pytest
import pandas as pd
import numpy as np
import logging

import topas_pipeline.metrics as metrics

# Mock logger to avoid unnecessary logging during tests
logger = logging.getLogger(__name__)


class TestComputeMetrics:
    def test_compute_metrics_with_existing_files_fp(self, mocker):
        # Mocking dependencies
        mocker.patch(
            "topas_pipeline.metrics.check_measures_computed", return_value=False
        )
        mocker.patch(
            "topas_pipeline.metrics.clinical_process.read_annotation_files",
            return_value=(
                pd.DataFrame(
                    {
                        "TOPAS_score": ["basket1", "basket2"],
                        "TOPAS_score_weights": [0.5, 1.0],
                        "TOPAS_subscore": ["sub_basket1", "sub_basket2"],
                        "TOPAS_subscore_weights": [0.25, 0.8],
                    }
                ),
                pd.DataFrame(),
            ),
        )
        measure_input = None
        measures_without_ref = None

        mock_save_measures = mocker.patch("topas_pipeline.metrics.save_measures")

        results_folder = Path("/path/to/results")
        debug = False
        data_types = ["fp"]

        metrics.compute_metrics(results_folder, debug, data_types)

        # Assertions to ensure the functions were called correctly
        assert metrics.check_measures_computed.call_count == len(data_types)
        assert metrics.clinical_process.read_annotation_files.call_count == len(
            data_types
        )
        assert metrics.save_measures.call_count == len(data_types)

        assert list(mock_save_measures.call_args[0][2].keys()) == [
            "rank",
            "fold_change",
            "z-score",
            "p-value",
            "occurrence",
        ]


class TestCheckMeasuresComputed:
    def test_all_measure_files_exist_fp(self, mocker):
        results_folder = "test_results"
        data_type = "fp"
        measure_names = [
            ("measure1", "measure1_short_name"),
            ("measure2", "measure2_short_name"),
        ]

        mocker.patch("os.path.exists", return_value=True)

        result = metrics.check_measures_computed(
            results_folder, data_type, measure_names
        )

        assert result is True

    def test_all_measure_files_exist_pp(self, mocker):
        results_folder = "test_results"
        data_type = "pp"
        measure_names = [
            ("measure1", "measure1_short_name"),
            ("measure2", "measure2_short_name"),
        ]

        mocker.patch("os.path.exists", return_value=True)

        result = metrics.check_measures_computed(
            results_folder, data_type, measure_names
        )

        assert result is True


class TestReadMeasures:
    # Correctly reads measures from existing TSV files
    def test_reads_measures_from_existing_tsv_files(self, mocker):
        mocker.patch("os.path.exists", return_value=True)
        mocker.patch(
            "pandas.read_csv",
            return_value=pd.DataFrame({"col1": [1, 2], "col2": [3, 4]}),
        )
        mocker.patch(
            "topas_pipeline.metrics.get_data_type_long", return_value="long_data_type"
        )
        mocker.patch(
            "topas_pipeline.metrics.utils.get_index_cols", return_value="index_col"
        )

        results_folder = "some/folder"
        data_type = "type"
        measure_names = [("m1", "measure1")]

        result = metrics.read_measures(results_folder, data_type, measure_names)

        expected_measure1_df = pd.DataFrame({"col1": [1, 2], "col2": [3, 4]})
        pd.testing.assert_frame_equal(result["measure1"], expected_measure1_df)

    def test_returns_empty_dict_one_file_missing(self, mocker):
        mocker.patch("os.path.exists", side_effect=[True, False])
        mocker.patch(
            "pandas.read_csv",
            return_value=pd.DataFrame({"col1": [1, 2], "col2": [3, 4]}),
        )
        mocker.patch(
            "topas_pipeline.metrics.get_data_type_long", return_value="long_data_type"
        )
        mocker.patch(
            "topas_pipeline.metrics.utils.get_index_cols", return_value="index_col"
        )

        results_folder = "some/folder"
        data_type = "type"
        measure_names = [("m1", "measure1"), ("m2", "measure1")]

        result = metrics.read_measures(results_folder, data_type, measure_names)

        assert result == {}


class TestSaveMeasures:
    def test_save_measures_correct_filenames(self, mocker):
        # Mocking the get_data_type_long function
        mocker.patch(
            "topas_pipeline.metrics.get_data_type_long", return_value="long_data_type"
        )

        # Mocking the to_csv method of pandas DataFrame
        mock_to_csv = mocker.patch.object(pd.DataFrame, "to_csv")

        results_folder = "test_folder"
        measure_names = [("measure1", "m1"), ("measure2", "m2")]
        measures = {
            "m1": pd.DataFrame({"col1": [1, 2], "col2": [3, 4]}),
            "m2": pd.DataFrame({"col1": [5, 6], "col2": [7, 8]}),
        }
        data_type = "type"

        metrics.save_measures(results_folder, measure_names, measures, data_type)

        expected_calls = [
            mocker.call(
                os.path.join(results_folder, "long_data_type_measures_measure1.tsv"),
                sep="\t",
                float_format="%.6g",
            ),
            mocker.call(
                os.path.join(results_folder, "long_data_type_measures_measure2.tsv"),
                sep="\t",
                float_format="%.6g",
            ),
        ]

        mock_to_csv.assert_has_calls(expected_calls, any_order=True)


class TestLooMedian:
    def test_loo_median_odd(self):
        input_matrix = np.array([[7.0, 6.0, 5.0, 4.1, 3.0, 2.0, 1.0]])

        expected_result = np.array([[3.55, 3.55, 3.55, 4.0, 4.55, 4.55, 4.55]])
        np.testing.assert_almost_equal(
            expected_result, metrics._loo_median(input_matrix)
        )

    def test_loo_median_odd_with_nan(self):
        input_matrix = np.array([[7.0, 6.0, 5.0, np.nan, 4.1, 3.0, 2.0, 1.0]])

        expected_result = np.array([[3.55, 3.55, 3.55, np.nan, 4.0, 4.55, 4.55, 4.55]])
        np.testing.assert_almost_equal(
            expected_result, metrics._loo_median(input_matrix)
        )

    def test_loo_median_even(self):
        input_matrix = np.array([[2.0, 7.0, 5.0, 3.0, 6.0, 1.0]])

        expected_result = np.array([[5.0, 3.0, 3.0, 5.0, 3.0, 5.0]])
        np.testing.assert_almost_equal(
            expected_result, metrics._loo_median(input_matrix)
        )

    def test_loo_median_two_values(self):
        input_matrix = np.array([[2.0, 7.0]])

        expected_result = np.array([[7.0, 2.0]])
        np.testing.assert_almost_equal(
            expected_result, metrics._loo_median(input_matrix)
        )

    def test_loo_median_one_value(self):
        input_matrix = np.array([[2.0]])

        with pytest.raises(IndexError):
            metrics._loo_median(input_matrix)

    def test_loo_median_one_value_and_nan(self):
        input_matrix = np.array([[2.0, np.nan]])

        expected_result = np.array([[np.nan, np.nan]])
        np.testing.assert_almost_equal(
            expected_result, metrics._loo_median(input_matrix)
        )


class TestLooStd:
    def test_loo_std_odd(self):
        input_matrix = np.array([[7.0, 6.0, 5.0, 4.1, 3.0, 2.0, 1.0]])

        expected_result = np.array(
            [
                [
                    1.8766104,
                    2.1637159,
                    2.3184046,
                    2.3664319,
                    2.3155273,
                    2.1575449,
                    1.8659225,
                ]
            ]
        )
        np.testing.assert_almost_equal(expected_result, metrics._loo_std(input_matrix))

    def test_loo_std_even(self):
        input_matrix = np.array([[2.0, 7.0, 5.0, 3.0, 6.0, 1.0]])

        expected_result = np.array(
            [[2.4083189, 2.0736441, 2.5884358, 2.5884358, 2.4083189, 2.0736441]]
        )
        np.testing.assert_almost_equal(expected_result, metrics._loo_std(input_matrix))

    def test_loo_std_even_with_nan(self):
        input_matrix = np.array([[2.0, 7.0, 5.0, np.nan, 3.0, 6.0, 1.0]])

        expected_result = np.array(
            [[2.4083189, 2.0736441, 2.5884358, np.nan, 2.5884358, 2.4083189, 2.0736441]]
        )
        np.testing.assert_almost_equal(expected_result, metrics._loo_std(input_matrix))

    def test_loo_std_two_values(self):
        input_matrix = np.array([[2.0, 7.0]])

        expected_result = np.array([[0.0, 0.0]])
        np.testing.assert_almost_equal(expected_result, metrics._loo_std(input_matrix))

    def test_loo_std_one_value(self):
        input_matrix = np.array([[2.0]])

        expected_result = np.array([[0.0]])
        np.testing.assert_almost_equal(expected_result, metrics._loo_std(input_matrix))

    def test_loo_std_one_value_and_nan(self):
        input_matrix = np.array([[2.0, np.nan]])

        expected_result = np.array([[0, np.nan]])
        np.testing.assert_almost_equal(expected_result, metrics._loo_std(input_matrix))


@pytest.fixture
def metrics_instance():
    # Sample DataFrame for testing
    data = {
        "protein1": [1.0, 2.0, 3.0, 4.0],
        "protein2": [2.0, 3.0, 4.0, np.nan],
        "protein3": [3.0, 4.0, 5.0, 6.0],
    }
    df = pd.DataFrame(data, index=["patient1", "patient2", "patient3", "patient4"])
    return metrics.Metrics(df.T)


def test_get_rank(metrics_instance):
    expected_rank = pd.DataFrame(
        {
            "rank_patient1": [4.0, 3.0, 4.0],
            "rank_patient2": [3.0, 2.0, 3.0],
            "rank_patient3": [2.0, 1.0, 2.0],
            "rank_patient4": [1.0, np.nan, 1.0],
            "rank_max": [4, 3, 4],
        },
        index=["protein1", "protein2", "protein3"],
    )
    result = metrics.Metrics.get_rank(metrics_instance.df)
    pd.testing.assert_frame_equal(result, expected_rank)


def test_get_fold_change(metrics_instance):
    # Define expected results manually or calculate based on the leave-one-out approach
    expected_fold_change = pd.DataFrame(
        {
            "fc_patient1": [0.01, 0.01818181818181818, 0.01],  # 0.01 = 10**(1.0-3.0)
            "fc_patient2": [0.1, 0.19801980198019803, 0.1],
            "fc_patient3": [10.0, 18.181818181818183, 10.0],
            "fc_patient4": [100.0, np.nan, 100.0],
        },
        index=["protein1", "protein2", "protein3"],
    )
    result = metrics.Metrics.get_fold_change(metrics_instance.df)
    pd.testing.assert_frame_equal(result, expected_fold_change)


def test_get_zscore(metrics_instance):
    # Manually calculate expected z-scores for verification
    expected_zscore = pd.DataFrame(
        {
            "zscore_patient1": [-2.0, -2.1213203435596424, -2.0],
            "zscore_patient2": [
                -0.6546536707079772,
                0.0,
                -0.6546536707079772,
            ],
            "zscore_patient3": [
                0.6546536707079772,
                2.1213203435596424,
                0.6546536707079772,
            ],
            "zscore_patient4": [2.0, np.nan, 2.0],
        },
        index=["protein1", "protein2", "protein3"],
    )
    result = metrics.Metrics.get_zscore(metrics_instance.df)
    pd.testing.assert_frame_equal(result, expected_zscore)


def test_get_pvalues(metrics_instance):
    z_scores = metrics.Metrics.get_zscore(metrics_instance.df)
    expected_pvalues = pd.DataFrame(
        {
            "pvalue_zscore_patient1": [
                0.04550026389635839,
                0.03389485352468927,
                0.04550026389635839,
            ],
            "pvalue_zscore_patient2": [
                0.5126907602619233,
                1.0,
                0.5126907602619233,
            ],
            "pvalue_zscore_patient3": [
                0.5126907602619233,
                0.03389485352468927,
                0.5126907602619233,
            ],
            "pvalue_zscore_patient4": [
                0.04550026389635839,
                np.nan,
                0.04550026389635839,
            ],
        },
        index=["protein1", "protein2", "protein3"],
    )
    result = metrics.Metrics.get_pvalues(z_scores)
    pd.testing.assert_frame_equal(result, expected_pvalues)


def test_calc(metrics_instance):
    metrics_instance.calc()
    assert "rank" in metrics_instance.metrics_df
    assert "fold_change" in metrics_instance.metrics_df
    assert "z-score" in metrics_instance.metrics_df
    assert "p-value" in metrics_instance.metrics_df
    assert "occurrence" in metrics_instance.metrics_df
    assert metrics_instance.metrics_df["occurrence"].shape == (3, 1)


if __name__ == "__main__":
    pytest.main()
