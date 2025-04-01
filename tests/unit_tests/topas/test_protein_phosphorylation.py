import pandas as pd

from topas_pipeline.topas import protein_phosphorylation


class TestProteinPhosphoScoring:
    def test_function_runs_with_valid_input(self, mocker):
        preprocessed_protein_df = pd.DataFrame(
            {
                "Gene names": ["GeneA", "GeneB"],
                "pat_col1": [1.0, 2.0],
                "pat_col2": [3.0, 4.0],
            }
        )

        mocker.patch("os.path.exists", return_value=False)
        mocker.patch("pandas.DataFrame.to_csv")
        mock_makedirs = mocker.patch("os.makedirs")
        protein_phosphorylation.protein_phospho_scoring(
            "dummy_results_folder", preprocessed_protein_df
        )

        mock_makedirs.assert_called_once_with("dummy_results_folder/protein_results")


class TestProteinScorePreprocess:
    # Returns preprocessed data if file already exists
    def test_filepath_does_not_exist_replica(self, mocker):
        mocker.patch("os.path.exists", return_value=False)
        mocker.patch(
            "pandas.read_csv",
            side_effect=[
                pd.DataFrame(
                    {
                        "zscore_patient1": [1.0],
                        "zscore_patient2": [2.0],
                        "Gene names": ["Gene1"],
                        "Modified sequence": ["Seq1"],
                        "Proteins": ["Prot1"],
                    }
                )
            ],
        )
        mocker.patch("pandas.DataFrame.to_csv")

        results_folder = "/path/to/results"
        result = protein_phosphorylation.protein_score_preprocess(results_folder)

        expected_result = pd.DataFrame(
            {
                "Gene names": ["Gene1"],
                "Proteins": ["Prot1"],
                "pat_patient1": [1.0],
                "pat_patient2": [2.0],
                "Peptide count": [2],
                "Peptide occurrence": ["2/2"],
            },
            index=pd.Series(["Seq1"], name="Modified sequence"),
        )

        pd.testing.assert_frame_equal(result, expected_result, check_dtype=False)
