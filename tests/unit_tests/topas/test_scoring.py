import os

import pytest
import pandas as pd
import numpy as np

from topas_pipeline.topas import scoring
import topas_pipeline.topas.scoring


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


class TestSumWeightedZScores:
    def test_filter_columns(self):
        data = {
            "kinase": ["A", "A", "B", "B", "C", "C"],
            "weight_patient_1": [1.2, 0.8, 1.5, 1.1, 1.3, 0.9],
            "weight_patient_2": [0.9, 1.1, 1.4, 1.0, 1.3, 1.2],
            "weighted_patient_1": [120, 140, 130, 110, 115, 135],
            "weighted_patient_2": [80, 85, 78, 90, 88, 84],
            "pat_patient_1": [110, 130, 125, 100, 110, 120],
            "pat_patient_2": [70, 75, 68, 85, 78, 76],
        }

        df = pd.DataFrame(data)
        result = topas_pipeline.topas.scoring.sum_weighted_z_scores(df, by="kinase")

        expected_result = pd.DataFrame(
            {
                "kinase": ["A", "B", "C"],
                "pat_patient_1": [260, 240, 250],
                "pat_patient_2": [165, 168, 172],
            }
        )

        pd.testing.assert_frame_equal(result, expected_result)


class TestCalculateModifiedSequenceWeights:
    def test_correctly_sums_weights_for_each_modified_sequence(self):
        data = {
            "Modified sequence": ["seq1", "seq1", "seq2", "seq2", "seq2"],
            "pat_1": [1, 1, 3, 3, 3],
            "pat_2": [5, 5, 7, 7, 7],
            "weight_1": [1, 2, 3, 4, 5],
            "weight_2": [5, 6, 7, 8, 9],
            "annotation_column": ["A", "A", "B", "B", "A"],
        }
        df = pd.DataFrame(data)
        result = topas_pipeline.topas.scoring.calculate_modified_sequence_weights(
            df, "annotation_column"
        )
        expected_data = {
            "annotation_column": ["A", "A", "B"],
            "Modified sequence": ["seq1", "seq2", "seq2"],
            "weight_1": [3, 5, 7],
            "weight_2": [11, 9, 15],
            "pat_1": [1, 3, 3],
            "pat_2": [5, 7, 7],
        }
        expected_df = pd.DataFrame(expected_data)
        pd.testing.assert_frame_equal(result, expected_df)


class TestCapZscoresAndWeights:
    def test_cap_zscores_and_weights(self):
        data = {
            "pat_1": [5, -5, 3, -3],
            "pat_2": [4.5, -4.5, 2, -2],
            "weight_1": [0.5, 1.5, 0.8, 1.2],
            "weight_2": [0.5, 1.5, 0.8, 1.2],
        }
        df = pd.DataFrame(data)
        result = topas_pipeline.topas.scoring.cap_zscores_and_weights(df)

        expected_result = pd.DataFrame(
            {
                "pat_1": [4.0, -4, 3, -3],
                "pat_2": [4.0, -4, 2, -2],
                "weight_1": [0.5, 1.0, 0.8, 1.0],
                "weight_2": [0.5, 1.0, 0.8, 1.0],
            }
        )
        pd.testing.assert_frame_equal(result, expected_result)


class TestSecondLevelZScoring:
    def test_properly_normalizes_patient_data(self):
        data = {
            "kinase": ["A", "B", "C"],
            "pat_patient_1": [4, 8, 6],
            "pat_patient_2": [7, 5, 9],
        }
        df = pd.DataFrame(data)
        result = topas_pipeline.topas.scoring.second_level_z_scoring(df, by_column="kinase")

        expected_result = pd.DataFrame(
            {
                "kinase": ["A", "B", "C"],
                "pat_patient_1": [
                    -0.7071067811865476,
                    0.7071067811865476,
                    -0.7071067811865476,
                ],
                "pat_patient_2": [
                    0.7071067811865476,
                    -0.7071067811865476,
                    0.7071067811865476,
                ],
                "mean": [5.5, 6.5, 7.5],
                "stdev": [2.1213203435596424, 2.1213203435596424, 2.1213203435596424],
            }
        )
        pd.testing.assert_frame_equal(result, expected_result)


class TestCalculatePeptideOccurrence:
    def test_calculate_peptide_occurrence_returns_modified_dataframe(self):
        data = {"pat_1": [np.nan, 1.0, 2.0], "pat_2": [3.0, np.nan, 4.0]}
        df = pd.DataFrame(data)

        result_df = topas_pipeline.topas.scoring.calculate_peptide_occurrence(df)
        expected_result_df = pd.DataFrame(
            {
                "pat_1": [np.nan, 1.0, 2.0],
                "pat_2": [3.0, np.nan, 4.0],
                "Peptide count": [1, 1, 2],
                "Peptide occurrence": ["1/2", "1/2", "2/2"],
            }
        )

        pd.testing.assert_frame_equal(result_df, expected_result_df, check_dtype=False)


class TestTopasScorePreprocess:
    def test_preprocess_data_when_file_does_not_exist(self, mocker):
        # Mocking os.path.exists to always return False
        mocker.patch("os.path.exists", return_value=False)

        # Mocking pd.read_csv to return a sample dataframe
        sample_patients_proteins = pd.DataFrame(
            {
                "Gene names": ["Gene1", "Gene2"],
                "Modified sequence": ["Seq1", "Seq2"],
                "Site positions": ["1", "2"],
                "Site positions identified (MQ)": ["1;2", "2"],
            }
        ).set_index('Modified sequence')
        sample_patients_zscores = pd.DataFrame(
            {
                "Gene names": ["Gene1", "Gene2"],
                "Modified sequence": ["Seq1", "Seq2"],
                "Site positions": ["1", "2"],
                "Site positions identified (MQ)": ["1;2", "2"],
                "pat_1": [0.1, 0.2],
                "pat_2": [0.3, 0.4],
            }
        )
        mocker.patch(
            "pandas.read_csv",
            side_effect=[sample_patients_proteins, sample_patients_zscores],
        )

        # Mocking to_csv to do nothing
        mocker.patch("pandas.DataFrame.to_csv")

        results_folder = "/path/to/results"
        result = topas_pipeline.topas.scoring.topas_score_preprocess(results_folder)

        expected_result = pd.DataFrame(
            {
                "Modified sequence": ["Seq1", "Seq1", "Seq2"],
                "Gene names": ["Gene1", "Gene1", "Gene2"],
                "pat_1": [0.1, 0.1, 0.2],
                "pat_2": [0.3, 0.3, 0.4],
                "All site positions": ["1", "1", "2"],
                "Site positions": ["1", "2", "2"],
                "Peptide count": [2, 2, 2],
                "Peptide occurrence": ["2/2", "2/2", "2/2"],
            }
        )

        pd.testing.assert_frame_equal(result, expected_result, check_dtype=False)


class TestCalculatePerPatientTargets:
    def test_group_and_count_patient_columns(self):
        data = {
            "kinase": ["A", "A", "B", "B"],
            "pat_1": [1, 2, 3, 4],
            "pat_2": [5, 6, 7, 8],
        }
        df = pd.DataFrame(data)
        result = topas_pipeline.topas.scoring.calculate_per_patient_targets(df, "kinase")

        expected_data = {"kinase": ["A", "B"], "targets_1": [2, 2], "targets_2": [2, 2]}
        expected_df = pd.DataFrame(expected_data)

        pd.testing.assert_frame_equal(result, expected_df)


class TestGetTargetSpace:
    # Correctly groups annotated peptides by the specified column
    def test_correct_grouping_by_column(self):
        data = {
            "kinase": ["A", "A", "B", "B"],
            "Modified sequence": ["seq1", "seq2", "seq1", "seq3"],
        }
        annotated_peptides_df = pd.DataFrame(data)
        scored_peptides_df = pd.DataFrame(
            {
                "kinase": ["A", "A", "B", "B"],
                "pat_1": [1, 2, 3, 4],
                "pat_2": [5, 6, 7, 8],
            }
        )

        result = topas_pipeline.topas.scoring.get_target_space(
            annotated_peptides_df, "kinase", scored_peptides_df
        )
        expected = pd.DataFrame(
            {
                "kinase": ["A", "B"],
                "No. of total targets": [2, 2],
                "targets_1": [2, 2],
                "targets_2": [2, 2],
            }
        )

        pd.testing.assert_frame_equal(result, expected)


def test_calculate_local_peptide_weight():
    # Sample input dataframe
    df = pd.DataFrame(
        {
            "Modified sequence": ["A", "B", "C", "D"],
            "pat_1": [10, 20, 30, 40],
            "pat_2": [5, 10, 15, np.nan],
            "Peptide count": [1, 2, 3, 4],
            "Site positions": ["X_S123", "Y_T234", "X_S123", "Y_T234"],
        }
    )

    # Expected output dataframe
    expected_df = pd.DataFrame(
        {
            "Modified sequence": ["A", "B", "C", "D"],
            "pat_1": [10, 20, 30, 40],
            "pat_2": [5, 10, 15, np.nan],
            "Peptide count": [1, 2, 3, 4],
            "Site positions": ["X_S123", "Y_T234", "X_S123", "Y_T234"],
            "weight_1": [1 / 4, 2 / 6, 3 / 4, 4 / 6],
            "weight_2": [1 / 4, 1, 3 / 4, np.nan],
        }
    )

    # Call the function
    result_df = topas_pipeline.topas.scoring.calculate_psite_weights(df)

    # Assert that the result matches the expected output
    pd.testing.assert_frame_equal(result_df, expected_df, check_like=True)


def test_calculate_weighted_zscores():
    df = pd.DataFrame(
        {
            "Modified sequence": ["A", "B", "C", "D"],
            "pat_1": [10, 20, 30, 40],
            "pat_2": [5, 10, 15, np.nan],
            "Peptide count": [1, 2, 3, 4],
            "Site positions": ["X_S123", "Y_T234", "X_S123", "Y_T234"],
            "weight_1": [1 / 4, 2 / 6, 3 / 4, 4 / 6],
            "weight_2": [1 / 4, 1, 3 / 4, np.nan],
        }
    )

    expected_df = pd.DataFrame(
        {
            "Modified sequence": ["A", "B", "C", "D"],
            "pat_1": [10, 20, 30, 40],
            "pat_2": [5, 10, 15, np.nan],
            "Peptide count": [1, 2, 3, 4],
            "Site positions": ["X_S123", "Y_T234", "X_S123", "Y_T234"],
            "weight_1": [1 / 4, 2 / 6, 3 / 4, 4 / 6],
            "weight_2": [1 / 4, 1, 3 / 4, np.nan],
            "weighted_1": [2.5, 20 / 3, 22.5, 80 / 3],
            "weighted_2": [5 / 4, 10, 45 / 4, np.nan],
        }
    )

    # Call the function
    result_df = topas_pipeline.topas.scoring.calculate_weighted_z_scores(df)

    # Assert that the result matches the expected output
    pd.testing.assert_frame_equal(result_df, expected_df, check_like=True)
