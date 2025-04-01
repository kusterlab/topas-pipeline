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


