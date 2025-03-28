import numpy as np
import pandas as pd
import pytest

from topas_pipeline.topas import topas


class TestComputeTopasScores:
    # Computes TOPAS scores correctly given valid input files and paths
    def test_compute_topas_scores_valid_inputs(self, topas_annotation_df, mocker):
        # Mocking dependencies
        mocker.patch(
            "topas_pipeline.sample_metadata.load",
            return_value=pd.DataFrame(
                {
                    "Sample name": ["Sample1", "Sample2", "Sample3"],
                    "Histologic subtype": ["subtype1", "subtype1", "subtype1"],
                }
            ),
        )
        mocker.patch(
            "topas_pipeline.topas.annotation.read_topas_annotations",
            return_value=topas_annotation_df,
        )
        mocker.patch(
            "topas_pipeline.topas.scoring.load_z_scores_fp",
            return_value=pd.DataFrame(
                {
                    "pat_Sample1": [0.5, 0.8],
                    "pat_Sample2": [0.6, 0.9],
                    "pat_Sample3": [0.7, 1.0],
                },
                index=pd.Series(["Gene1", "Gene2"], name="Gene names"),
            ),
        )
        mocker.patch(
            "topas_pipeline.topas.scoring.load_z_scores_pp",
            return_value=pd.DataFrame(
                {
                    "pat_Sample1": [1.0, 2.0],
                    "pat_Sample2": [3, 4],
                    "pat_Sample3": [5, 6],
                },
                index=pd.Series(["seq1", "seq2"], name="Modified sequence"),
            ),
        )
        mocker.patch(
            "topas_pipeline.topas.scoring.load_protein_phosphorylation",
            return_value=pd.DataFrame(
                {
                    "pat_Sample1": [0.5, 0.8],
                    "pat_Sample2": [0.6, 0.9],
                    "pat_Sample3": [0.7, 1.0],
                },
                index=pd.Series(["Gene1", "Gene2"], name="Gene names"),
            ),
        )
        mocker.patch(
            "topas_pipeline.topas.scoring.load_kinase_scores",
            return_value=pd.DataFrame(
                {
                    "pat_Sample1": [0.5, 0.8],
                    "pat_Sample2": [0.6, 0.9],
                    "pat_Sample3": [0.7, 1.0],
                },
                index=pd.Series(["Kinase1", "Kinase2"], name="PSP Kinases"),
            ),
        )
        mock_save_topas_scores = mocker.patch(
            "topas_pipeline.topas.topas.save_topas_scores"
        )
        mocker.patch("topas_pipeline.topas.topas.save_rtk_scores_w_metadata")
        mocker.patch("pandas.DataFrame.to_csv")

        # Define inputs
        results_folder = "results"
        metadata_file = "metadata.csv"
        topas_file = "topas_scores.csv"
        topas_results_folder = "topas_results"

        # Call the function
        topas.compute_topas_scores(
            results_folder=results_folder,
            metadata_file=metadata_file,
            topas_annotation_file=topas_file,
            topas_results_folder=topas_results_folder,
        )

        # first [0]: first call
        # second [0]: first of (args, kwargs) tuple
        # third [0]: first argument
        result = mock_save_topas_scores.call_args_list[0][0][0]

        expected_result = pd.DataFrame(
            {"topas1": [0.5, 0.6, 0.7]},
            index=pd.Series(
                ["pat_Sample1", "pat_Sample2", "pat_Sample3"], name="index"
            ),
        )

        pd.testing.assert_frame_equal(result, expected_result)

        result_zscored = mock_save_topas_scores.call_args_list[1][0][0]

        expected_result_zscored = pd.DataFrame(
            {"topas1": [-2.121320343559642, 0.0, 2.121320343559642]},
            index=pd.Series(
                ["pat_Sample1", "pat_Sample2", "pat_Sample3"], name="index"
            ),
        )

        pd.testing.assert_frame_equal(result_zscored, expected_result_zscored)


@pytest.fixture
def topas_annotation_df():
    return pd.DataFrame(
        {
            "TOPAS_score": ["topas1", "topas2"],
            "TOPAS_subscore": ["subtopas1", "subtopas2"],
            "Gene names": ["Gene1", "Gene2"],
            "weight": [1, np.nan],
            "group": ["group1", "OTHER"],
            "Scoring rule": ["highest z-score", "highest z-score"],
            "TOPAS_subscore_level": ["level1", "level1"],
        }
    )


class TestReadTopasScores:
    def test_reads_topas_scores_from_valid_file_path(self, mocker):
        # Mocking os.path.exists to return True
        mocker.patch("os.path.exists", return_value=True)
        # Mocking pd.read_csv to return a sample DataFrame
        sample_df = pd.DataFrame({"Sample": ["A", "B"], "Score": [1, 2]}).set_index(
            "Sample"
        )
        mocker.patch("pandas.read_csv", return_value=sample_df)

        results_folder = "some/folder"
        z_scored = False

        result = topas.read_topas_scores(results_folder, z_scored)

        expected_result = pd.DataFrame(
            {"Score": [1, 2]}, index=pd.Series(["A", "B"], name="Sample")
        ).T

        pd.testing.assert_frame_equal(result, expected_result)


class TestSaveRtkScoresWMetadata:
    def test_correct_subsetting_to_rtk(self, mocker):
        TOPAS_CATEGORIES = {"RTK1": "value1", "RTK2": "value2"}
        mocker.patch("topas_pipeline.topas.topas.TOPAS_CATEGORIES", TOPAS_CATEGORIES)
        mock_to_csv = mocker.patch("pandas.DataFrame.to_csv", autospec=True)

        # Sample data
        topas_scores_data = {
            "Sample name": ["Sample1", "Sample2"],
            "RTK1": [0.1, 0.4],
            "RTK2": [0.2, 0.5],
            "Other": [0.3, 0.6],
        }
        metadata_data = {
            "Sample name": ["Sample1", "Sample2"],
            "code_oncotree": ["M1", "M3"],
            "Meta2": ["M2", "M4"],
        }

        topas_scores_df = pd.DataFrame(topas_scores_data).set_index("Sample name")
        metadata_df = pd.DataFrame(metadata_data)

        # Output file path
        out_file = "output.tsv"

        # Call the function
        topas.save_rtk_scores_w_metadata(topas_scores_df, metadata_df, out_file)

        result = mock_to_csv.call_args[0][0]
        expected_result = pd.DataFrame(
            {"RTK1": [0.1, 0.4], "RTK2": [0.2, 0.5], "code_oncotree": ["M1", "M3"]},
            index=pd.Series(["Sample1", "Sample2"], name="Sample name"),
        )
        print(result.to_dict(orient="list"))
        pd.testing.assert_frame_equal(result, expected_result)


class TestSaveTopasScores:
    def test_renames_index_by_removing_score_prefix(self, mocker):
        data = {"A": [1, 2], "B": [3, 4]}
        index = ["score_1", "score_2"]
        df = pd.DataFrame(data, index=index)

        mock_to_csv = mocker.patch("pandas.DataFrame.to_csv", autospec=True)

        out_file = "out_file.tsv"
        topas.save_topas_scores(df, out_file)

        result = mock_to_csv.call_args[0][0]

        expected_result = pd.DataFrame(
            {"A": [1, 2], "B": [3, 4]}, index=pd.Index(["1", "2"], name="Sample")
        )

        pd.testing.assert_frame_equal(result, expected_result)


class TestSaveTopasScoresLongFormat:
    def test_converts_wide_to_long_format_correctly(self, mocker):
        data = {
            "Sample": ["S1", "S2"],
            "FP.RTK Scores": [0.1, 0.2],
            "PP.Main Scores": [0.3, 0.4],
        }
        topas_scores_df = pd.DataFrame(data)
        out_file = "out_file.tsv"

        mock_to_csv = mocker.patch("pandas.DataFrame.to_csv", autospec=True)

        topas.save_topas_scores_long_format(topas_scores_df, out_file)

        result = mock_to_csv.call_args[0][0]

        expected_data = {
            "Sample": ["S1", "S2", "S1", "S2"],
            "variable": ["RTK Scores", "RTK Scores", "Main Scores", "Main Scores"],
            "value": [0.1, 0.2, 0.3, 0.4],
            "Data type": ["FP", "FP", "PP", "PP"],
            "Basket type": ["RTK", "RTK", "Main", "Main"],
        }
        expected_df = pd.DataFrame(expected_data)

        pd.testing.assert_frame_equal(result, expected_df)
