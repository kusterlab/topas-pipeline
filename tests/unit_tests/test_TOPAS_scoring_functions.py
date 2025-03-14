import pandas as pd
import numpy as np

import topas_pipeline.TOPAS_scoring_functions as scoring_functions


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
        result = scoring_functions.calculate_modified_sequence_weights(
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
        result = scoring_functions.cap_zscores_and_weights(df)

        expected_result = pd.DataFrame(
            {
                "pat_1": [4.0, -4, 3, -3],
                "pat_2": [4.0, -4, 2, -2],
                "weight_1": [0.5, 1.0, 0.8, 1.0],
                "weight_2": [0.5, 1.0, 0.8, 1.0],
            }
        )
        pd.testing.assert_frame_equal(result, expected_result)


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
        result = scoring_functions.sum_weighted_z_scores(df, by="kinase")

        expected_result = pd.DataFrame(
            {
                "kinase": ["A", "B", "C"],
                "pat_patient_1": [260, 240, 250],
                "pat_patient_2": [165, 168, 172],
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
        result = scoring_functions.second_level_z_scoring(df, by_column="kinase")

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

        result_df = scoring_functions.calculate_peptide_occurrence(df)
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
                "Other": ["Data1", "Data2"],
            }
        )
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
        result = scoring_functions.topas_score_preprocess(results_folder)

        expected_result = pd.DataFrame(
            {
                "Modified sequence": ["Seq1", "Seq1", "Seq2"],
                "Gene names": ["Gene1", "Gene1", "Gene2"],
                "pat_1": [0.1, 0.1, 0.2],
                "pat_2": [0.3, 0.3, 0.4],
                "All site positions": ["1", "1", "2"],
                "Site positions": ["1", "2", "2"],
                "Other": ["Data1", "Data1", "Data2"],
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
        result = scoring_functions.calculate_per_patient_targets(df, "kinase")

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

        result = scoring_functions.get_target_space(
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
    result_df = scoring_functions.calculate_psite_weights(df)

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
    result_df = scoring_functions.calculate_weighted_z_scores(df)

    # Assert that the result matches the expected output
    pd.testing.assert_frame_equal(result_df, expected_df, check_like=True)
