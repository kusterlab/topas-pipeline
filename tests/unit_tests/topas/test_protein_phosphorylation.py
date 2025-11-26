import pandas as pd

from topas_pipeline.topas import protein_phosphorylation


class TestProteinPhosphoScoring:
    def test_function_runs_with_valid_input(self, mocker):
        # Mock Path.is_file to pretend the output file does NOT exist
        mocker.patch(
            "pathlib.Path.is_file",
            return_value=False,
        )

        # Mock get_cohort_intensities_df to return a minimal dataframe
        mock_cohort_df = pd.DataFrame(
            {
                "Gene names": ["GeneA", "GeneB"],
                "pat_col1": [1.0, 2.0],
                "pat_col2": [3.0, 4.0],
            }
        ).set_index("Gene names")

        mocker.patch(
            "topas_pipeline.preprocess.phospho_grouping.read_cohort_intensities_df",
            return_value=mock_cohort_df,
        )

        # Mock compute_substrate_phosphorylation_scores
        mock_compute = mocker.patch(
            "topas_pipeline.topas.ck_substrate_phosphorylation.compute_substrate_phosphorylation_scores",
            return_value="dummy_scores",
        )

        # Mock save_scores
        mock_save = mocker.patch(
            "topas_pipeline.topas.ck_substrate_phosphorylation.save_scores"
        )

        # Run function
        protein_phosphorylation.protein_phospho_scoring(
            "dummy_results_folder",
            metadata_file="dummy_meta.tsv",
        )

        # Verify scoring function called
        mock_compute.assert_called_once()

        # Verify scores saved twice (once normal, once metadata)
        assert mock_save.call_count == 2
