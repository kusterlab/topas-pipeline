import os

import pytest
import pandas as pd
import numpy as np

from topas_pipeline.topas import scoring


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

        result = scoring.read_topas_scores(results_folder, z_scored)

        expected_result = pd.DataFrame(
            {"Score": [1, 2]}, index=pd.Series(["A", "B"], name="Sample")
        ).T

        pd.testing.assert_frame_equal(result, expected_result)


class TestReadTopasSubscores:
    def test_reads_topas_subscore_files_correctly(self, mocker):
        # Mocking the get_paths_to_sub_basket_files function
        mocker.patch(
            "topas_pipeline.topas.scoring.get_paths_to_topas_subscore_files",
            return_value=["file1.tsv", "file2.tsv"],
        )

        # Mocking os.path.exists to always return True
        mocker.patch("os.path.exists", return_value=True)

        # Creating a sample dataframe to return when pd.read_csv is called
        sample_df = pd.DataFrame(
            {
                "index": [1, 2],
                "total_basket_score": [0.5, 0.6],
                "score1": [0.1, 0.2],
                "score2": [0.3, 0.4],
            }
        ).set_index("index")

        # Mocking pd.read_csv to return the sample dataframe
        mocker.patch("pandas.read_csv", return_value=sample_df)

        results_folder = "dummy_folder"
        result = scoring.read_topas_subscores(results_folder)

        # Expected dataframe after processing
        expected_df = pd.DataFrame(
            {1: [0.1, 0.3, 0.1, 0.3], 2: [0.2, 0.4, 0.2, 0.4]},
            index=pd.Series(["score1", "score2", "score1", "score2"], name=None),
        )
        expected_df.columns.name = "index"

        pd.testing.assert_frame_equal(result, expected_df)


class TestLoadZScoresFp:
    # Correctly reads z-scores from the specified results folder
    def test_correctly_reads_z_scores(self, mocker):
        # Mocking the dependencies
        measure1_df = pd.DataFrame(
            {"zscore_pat_patient1": [1, np.nan], "zscore_pat_patient2": [3, 4]}
        )
        mocker.patch(
            "topas_pipeline.metrics.read_measures",
            return_value={"z-score": measure1_df},
        )
        mocker.patch(
            "topas_pipeline.metrics.clinical_annotation.read_annotated_expression_file",
            return_value=(
                pd.DataFrame(
                    {
                        "Identification metadata patient1": ["", "detected in batch;"],
                        "Identification metadata patient2": ["", ""],
                    }
                )
            ),
        )

        results_folder = "some/results/folder"
        result = scoring.load_z_scores_fp(results_folder)

        expected_result = pd.DataFrame(
            {"pat_patient1": [1.0, -4], "pat_patient2": [3, 4]}
        )

        pd.testing.assert_frame_equal(result, expected_result)


class TestLoadZScoresPp:
    # Correctly reads z-scores from the specified results folder
    def test_correctly_reads_z_scores(self, mocker):
        # Mocking the dependencies
        measure1_df = pd.DataFrame(
            {
                "Modified sequence": ["seq1", "seq2"],
                "zscore_pat_patient1": [1, np.nan],
                "zscore_pat_patient2": [3, 4],
            }
        )
        mocker.patch(
            "topas_pipeline.metrics.read_measures",
            return_value={"z-score": measure1_df},
        )

        results_folder = "some/results/folder"
        result = scoring.load_z_scores_pp(results_folder)

        expected_result = pd.DataFrame(
            {"pat_patient1": [1.0, np.nan], "pat_patient2": [3, 4]},
            index=pd.Series(["seq1", "seq2"], name="Modified sequence"),
        )

        pd.testing.assert_frame_equal(result, expected_result)


class TestLoadProteinPhosphorylation:
    def test_loads_and_parses_tsv_correctly(self, mocker):
        mock_read_csv = mocker.patch("pandas.read_csv")
        mock_read_csv.return_value = pd.DataFrame(
            {
                "Gene names": ["Gene1", "Gene2", "Gene2;Gene3"],
                "pat_Sample1": [0.5, 0.8, 1.1],
                "pat_Sample2": [0.6, 0.9, 1.2],
            }
        ).set_index("Gene names")

        results_folder = "test_folder"
        result = scoring.load_protein_phosphorylation(results_folder)

        # Assertions
        mock_read_csv.assert_called_once_with(
            os.path.join(results_folder, "protein_results/protein_scores.tsv"),
            sep="\t",
            index_col="Gene names",
        )

        expected_result = pd.DataFrame(
            {"pat_Sample1": [0.5, 0.8], "pat_Sample2": [0.6, 0.9]},
            index=pd.Series(["Gene1", "Gene2"], name="Gene names"),
        )
        pd.testing.assert_frame_equal(result, expected_result)


class TestLoadKinaseScores:
    def test_loads_and_parses_tsv_correctly(self, mocker):
        # Mocking the pd.read_csv function
        mock_read_csv = mocker.patch("pandas.read_csv")
        mock_read_csv.return_value = pd.DataFrame(
            {
                "PSP Kinases": ["Kinase1", "Kinase2"],
                "pat_Sample1": [0.5, 0.8],
                "pat_Sample2": [0.6, 0.9],
            }
        ).set_index("PSP Kinases")

        results_folder = "test_folder"
        result = scoring.load_kinase_scores(results_folder)

        # Assertions
        mock_read_csv.assert_called_once_with(
            os.path.join(results_folder, "kinase_results/kinase_scores.tsv"),
            sep="\t",
            index_col="PSP Kinases",
        )

        expected_result = pd.DataFrame(
            {"pat_Sample1": [0.5, 0.8], "pat_Sample2": [0.6, 0.9]},
            index=pd.Series(["Kinase1", "Kinase2"], name="PSP Kinases"),
        )
        pd.testing.assert_frame_equal(result, expected_result)


class TestGetNumberIdentAnnotPerSample:
    def test_reads_csv_files_correctly(self, mocker):
        # Mocking pd.read_csv to return a sample annot_df DataFrame
        mock_df = pd.DataFrame(
            {
                "Gene names": ["gene1", "gene2"],
                "pat_1": [1, 2],
                "pat_2": [3, 4],
                "TOPAS_score": ["A", "B"],
            }
        )
        mocker.patch("pandas.read_csv", return_value=mock_df)

        results_folder = "test_folder"
        data_types = ["fp"]

        result = scoring.get_number_ident_annot_per_sample(
            results_folder, data_types
        )
        expected_result = pd.DataFrame(
            {"fp.num_identified": [2, 2], "fp.num_annotated": [2, 2]},
            index=pd.Series(["pat_1", "pat_2"]),
        )
        pd.testing.assert_frame_equal(result, expected_result)


class TestComputeTopasScores:
    # Computes TOPAS scores correctly given valid input files and paths
    def test_compute_topas_scores_valid_inputs(self, topas_annotation_df, mocker):
        # Mocking dependencies
        mocker.patch(
            "topas_pipeline.topas.scoring.sample_metadata.load",
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
        mocker.patch(
            "topas_pipeline.topas.scoring.get_number_ident_annot_per_sample",
            return_value=pd.DataFrame(
                {"fp.num_identified": [2, 2], "fp.num_annotated": [2, 2]},
                index=pd.Series(["pat_Sample1", "pat_Sample2"]),
            ),
        )
        mock_save_topas_scores = mocker.patch(
            "topas_pipeline.topas.scoring.save_topas_scores"
        )
        mocker.patch("topas_pipeline.topas.scoring.save_rtk_scores_w_metadata")
        mocker.patch("pandas.DataFrame.to_csv")

        # Define inputs
        results_folder = "results"
        metadata_file = "metadata.csv"
        topas_file = "topas_scores.csv"
        topas_results_folder = "topas_results"

        # Call the function
        scoring.compute_TOPAS_scores(
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


class TestSaveRtkScoresWMetadata:
    def test_correct_subsetting_to_rtk(self, mocker):
        TOPAS_CATEGORIES = {"RTK1": "value1", "RTK2": "value2"}
        mocker.patch("topas_pipeline.topas.scoring.TOPAS_CATEGORIES", TOPAS_CATEGORIES)
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
        scoring.save_rtk_scores_w_metadata(topas_scores_df, metadata_df, out_file)

        result = mock_to_csv.call_args[0][0]
        expected_result = pd.DataFrame(
            {"RTK1": [0.1, 0.4], "RTK2": [0.2, 0.5], "code_oncotree": ["M1", "M3"]},
            index=pd.Series(["Sample1", "Sample2"], name="Sample name"),
        )
        print(result.to_dict(orient="list"))
        pd.testing.assert_frame_equal(result, expected_result)


class TestCountSignificantTopasScores:
    def test_default_threshold(self):
        data = {
            "topas1": [2.1, 1.9, 2.5],
            "topas2": [1.8, 2.2, 2.0],
            "topas3": [2.0, 2.1, 1.7],
        }
        df_zscored = pd.DataFrame(data)

        result = scoring.count_significant_topas_scores(df_zscored)

        expected_data = {"num_significant_baskets": [2, 2, 2]}
        expected_df = pd.DataFrame(expected_data)

        pd.testing.assert_frame_equal(result, expected_df)

    def test_multiple_threshold(self):
        data = {
            "topas1": [2.1, 1.9, 2.5],
            "topas2": [1.8, 2.2, 2.0],
            "topas3": [2.0, 2.1, 1.6],
        }
        df_zscored = pd.DataFrame(data)

        result = scoring.count_significant_topas_scores(
            df_zscored, thresholds=[1.7, 2.0]
        )

        expected_data = {
            "num_significant_baskets": [3, 3, 2],
            "num_significant_baskets_z>=2.0": [2, 2, 2],
        }
        expected_df = pd.DataFrame(expected_data)

        pd.testing.assert_frame_equal(result, expected_df)


class TestGetSummedZscore:
    # Correctly calculates summed z-scores for non-ligand cases
    def test_summed_zscore_non_ligand(self, mocker):
        z_score_data = {"sample1": [1.5, -3.2, 4.5], "sample2": [-2.5, 3.0, -4.1]}
        topas_subscore_data = {
            "GENES": ["Gene4", "Gene1", "Gene2"],
            "weight": [0.5, 1.0, 2.0],
        }
        mocker.patch("pandas.DataFrame.to_csv")

        z_score_df = pd.DataFrame(
            z_score_data,
            index=pd.Series(["Gene1", "Gene2;Gene3", "Gene4"], name="Gene names"),
        )
        topas_subscore_df = pd.DataFrame(topas_subscore_data)

        result = scoring.get_summed_zscore(
            z_score_df, topas_subscore_df, "Gene names", "GENES", ligand=False
        )
        # 1.5*1 - 4.0 + 2.25*0.5 = -0.25
        # -2.5*1 + 4.0 + -2.05*0.5 = -0.55
        expected_result = pd.Series({"sample1": -0.25, "sample2": -0.55})

        pd.testing.assert_series_equal(result, expected_result)


class TestSaveTopasScores:
    def test_renames_index_by_removing_score_prefix(self, mocker):
        data = {"A": [1, 2], "B": [3, 4]}
        index = ["score_1", "score_2"]
        df = pd.DataFrame(data, index=index)

        mock_to_csv = mocker.patch("pandas.DataFrame.to_csv", autospec=True)

        out_file = "out_file.tsv"
        scoring.save_topas_scores(df, out_file)

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

        scoring.save_topas_scores_long_format(topas_scores_df, out_file)

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


class TestExtractTopasMemberZScores:
    def test_extracts_z_scores_correctly(self, topas_annotation_df, mocker):
        # Mocking the dependencies
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
        mocker.patch("os.makedirs")
        mock_to_csv = mocker.patch("pandas.DataFrame.to_csv", autospec=True)

        # Call the function
        scoring.extract_topas_member_z_scores("results_folder", "topas_file")

        result_gene1_level1 = mock_to_csv.call_args_list[0][0][0]
        expected_result_gene1_level1 = pd.DataFrame(
            {"Gene1 - level1": [0.5, 0.6, 0.7]},
            index=pd.Index(["pat_Sample1", "pat_Sample2", "pat_Sample3"]),
        )
        expected_result_gene1_level1.columns.name = "Gene names"
        pd.testing.assert_frame_equal(result_gene1_level1, expected_result_gene1_level1)

        result_gene2_level1 = mock_to_csv.call_args_list[1][0][0]
        expected_result_gene2_level1 = pd.DataFrame(
            {"Gene2 - level1": [0.8, 0.9, 1.0]},
            index=pd.Index(["pat_Sample1", "pat_Sample2", "pat_Sample3"]),
        )
        expected_result_gene2_level1.columns.name = "Gene names"
        pd.testing.assert_frame_equal(result_gene2_level1, expected_result_gene2_level1)


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
