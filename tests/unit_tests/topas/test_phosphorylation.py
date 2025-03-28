import pandas as pd

from topas_pipeline.topas import phosphorylation


class TestReadTopasSubscores:
    def test_reads_topas_subscore_files_correctly(self, mocker):
        # Mocking the get_paths_to_sub_basket_files function
        mocker.patch(
            "topas_pipeline.topas.phosphorylation.get_paths_to_topas_subscore_files",
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
        result = phosphorylation.read_topas_subscores(results_folder)

        # Expected dataframe after processing
        expected_df = pd.DataFrame(
            {1: [0.1, 0.3, 0.1, 0.3], 2: [0.2, 0.4, 0.2, 0.4]},
            index=pd.Series(["score1", "score2", "score1", "score2"], name=None),
        )
        expected_df.columns.name = "index"

        pd.testing.assert_frame_equal(result, expected_df)
