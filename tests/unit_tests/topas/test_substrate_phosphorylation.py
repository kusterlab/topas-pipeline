import pandas as pd
from topas_pipeline.topas.substrate_phosphorylation import kinase_scoring


class TestKinaseScoring:
    def test_kinase_scoring_valid_inputs(self, mocker):
        # Mocking the logger
        mock_logger = mocker.patch("topas_pipeline.topas.substrate_phosphorylation.logger.info")

        # Mocking the scoring functions
        mocker.patch(
            "topas_pipeline.topas.scoring.calculate_psite_weights",
            return_value=pd.DataFrame(columns=["PSP Kinases"]),
        )
        mocker.patch(
            "topas_pipeline.topas.scoring.calculate_modified_sequence_weights",
            return_value=pd.DataFrame(),
        )
        mocker.patch(
            "topas_pipeline.topas.scoring.cap_zscores_and_weights",
            return_value=pd.DataFrame(),
        )
        mocker.patch(
            "topas_pipeline.topas.scoring.calculate_weighted_z_scores",
            return_value=pd.DataFrame(),
        )
        mocker.patch(
            "topas_pipeline.topas.scoring.sum_weighted_z_scores",
            return_value=pd.DataFrame(),
        )
        mocker.patch(
            "topas_pipeline.topas.scoring.second_level_z_scoring",
            return_value=pd.DataFrame(columns=["PSP Kinases"]),
        )
        mocker.patch(
            "topas_pipeline.topas.scoring.get_target_space",
            return_value=pd.DataFrame(columns=["PSP Kinases", "No. of total targets"]),
        )

        # Mocking os.path.exists and os.makedirs
        mocker.patch("os.path.exists", return_value=False)
        mocker.patch("os.makedirs")

        # Mocking pandas to_csv method
        mocker.patch("pandas.DataFrame.to_csv")

        # Create a sample DataFrame for testing
        preprocessed_df = pd.DataFrame()

        # Call the function with valid inputs
        kinase_scoring("test_output_folder", preprocessed_df)

        # Assertions to ensure the function runs without errors
        assert mock_logger.called
