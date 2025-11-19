import io
from pathlib import Path

import pytest
import numpy as np
import pandas as pd

import topas_pipeline.preprocess.preprocess_tools as prep
from topas_pipeline.preprocess import sample_mapping
from topas_pipeline.sample_annotation import get_unique_batches
from topas_pipeline.identification_metadata import mark_num_peptides
from topas_pipeline.data_loaders import data_loader, tmt_loader


class TestCheckAnnot:
    # Loads sample annotation file correctly
    def test_loads_sample_annotation_file_correctly(self, mocker):
        mock_sample_annotation = mocker.patch(
            "topas_pipeline.sample_annotation.load_sample_annotation"
        )
        mock_filter_sample_annotation = mocker.patch(
            "topas_pipeline.sample_annotation.filter_sample_annotation"
        )
        mock_sample_annotation.return_value = mocker.Mock()
        mock_filter_sample_annotation.return_value = pd.DataFrame(
            columns=["Sample name", "Cohort", "Batch Name", "TMT Channel"]
        )
        mocker.patch(
            "pandas.read_excel",
            return_value=pd.DataFrame(columns=["Sample name"]),
        )

        sample_annotation_file = "sample_file.csv"
        metadata_annotation_file = "metadata_file.csv"
        in_metadata = mocker.Mock()

        result = prep.check_annot(
            sample_annotation_file, metadata_annotation_file, in_metadata
        )

        mock_sample_annotation.assert_called_once_with(sample_annotation_file)
        mock_filter_sample_annotation.assert_called_once()
        assert result is not None

        # Check for duplicated sample names

    def test_check_duplicated_sample_names(self, mocker):
        mock_load_sample_annotation = mocker.patch(
            "topas_pipeline.sample_annotation.load_sample_annotation"
        )
        mock_filter_sample_annotation = mocker.patch(
            "topas_pipeline.sample_annotation.filter_sample_annotation"
        )
        mock_logger_info = mocker.patch("topas_pipeline.preprocess_tools.logger.info")

        sample_annotation_file = "sample_file.csv"
        metadata_annotation_file = "metadata_file.csv"
        in_metadata = mocker.Mock()
        mocker.patch(
            "pandas.read_excel",
            return_value=pd.DataFrame(columns=["Sample name"]),
        )

        sample_annot_df_filtered = pd.DataFrame(
            {
                "Sample name": ["sample1", "sample2", "sample1"],
                "Batch Name": ["batch1", "batch1", "batch1"],
                "TMT Channel": [1, 2, 3],
            }
        )
        mock_load_sample_annotation.return_value = sample_annot_df_filtered
        mock_filter_sample_annotation.return_value = sample_annot_df_filtered

        with pytest.raises(
            ValueError,
            match=r"Duplicated sample\(s\) in sample annotation:",
        ):
            prep.check_annot(
                sample_annotation_file, metadata_annotation_file, in_metadata
            )

        mock_logger_info.assert_called_once()


class TestInMetadata:
    # Returns the sample string if it does not end with '-R2'
    def test_returns_sample_if_no_r2(self):
        sample = "CHD-123"
        result = prep.in_metadata(sample)
        assert result == sample

    # Sample string with multiple '-R2' occurrences
    def test_sample_with_multiple_r2(self):
        sample = "CHD-R2-R2"
        result = prep.in_metadata(sample)
        assert result == "CHD-R2"


class TestRemoveEmptyRefBatch:
    # Correctly removes 'Reporter intensity corrected' columns for empty batches based on sample annotation data
    def test_drop_reporter_intensity_corrected_columns_replica(self):
        data = {
            "Reporter intensity corrected 1 Batch3": [1, 2, 3],
            "Reporter intensity corrected 2 Batch210": [4, 5, 6],
            "Other column Batch5": [7, 8, 9],
        }
        df_with_ref = pd.DataFrame(data)
        sample_annotation_data = {
            "Batch Name": [3, 210],
            "Cohort": ["Cohort1", "Cohort 1"],
            "QC": ["passed", "failed"],
        }
        sample_annotation_df = pd.DataFrame(sample_annotation_data)

        result_df = prep.remove_ref_empty_batch(df_with_ref, sample_annotation_df)

        assert "Reporter intensity corrected 1 Batch3" in result_df.columns
        assert "Reporter intensity corrected 2 Batch210" not in result_df.columns
        assert "Other column Batch5" in result_df.columns


class TestGetFilesByType:
    # Correctly retrieves file paths based on provided data type and file type
    def test_correct_file_retrieval(self, mocker):
        mock_sample_annotation = mocker.patch("topas_pipeline.preprocess_tools.sample_annotation")
        mock_get_data_location = mocker.patch("topas_pipeline.preprocess_tools.get_data_location")
        mock_filter_evidence_files = mocker.patch(
            "topas_pipeline.preprocess_tools.filter_evidence_files"
        )

        sample_annotation_df = pd.DataFrame({"batch": ["batch1", "batch2"]})
        raw_data_location = "/path/to/raw/data"
        data_type = "typeA"
        file_type = "fileA"

        mock_sample_annotation.get_unique_batches.return_value = ["batch1", "batch2"]
        mock_get_data_location.return_value = [
            "/path/to/raw/data/fileA1",
            "/path/to/raw/data/fileA2",
        ]
        mock_filter_evidence_files.return_value = ["/path/to/raw/data/fileA1"]

        result = prep.get_files_by_type(
            raw_data_location, data_type, sample_annotation_df, file_type
        )

        assert result == ["/path/to/raw/data/fileA1"]
        mock_sample_annotation.get_unique_batches.assert_called_once_with(
            sample_annotation_df
        )
        mock_get_data_location.assert_called_once_with(
            raw_data_location, data_type, file_type=file_type
        )
        mock_filter_evidence_files.assert_called_once_with(
            ["/path/to/raw/data/fileA1", "/path/to/raw/data/fileA2"],
            data_type.upper(),
            ["batch1", "batch2"],
        )


class TestLogTransfromIntensities:
    # Log transformation is applied correctly to all intensity columns
    def test_log_transform_applied_correctly(self, mocker):
        data = {
            "Reporter intensity corrected 1": [100.0, 200.0, 300.0],
            "Reporter intensity corrected 2": [400.0, 500.0, 600.0],
            "Metadata column": [400, 500, 600],
        }
        df = pd.DataFrame(data)

        # Applying the log transformation
        result_df = prep.log_transform_intensities(df)

        # Expected DataFrame after log transformation
        expected_data = {
            "Reporter intensity corrected 1": np.log10([100.0, 200.0, 300.0]),
            "Reporter intensity corrected 2": np.log10([400.0, 500.0, 600.0]),
            "Metadata column": [400, 500, 600],
        }
        expected_df = pd.DataFrame(expected_data)

        # Asserting the transformation is correct
        pd.testing.assert_frame_equal(result_df, expected_df)


class TestSumPeptideIntensities:
    def test_aggregate_intensities_lfq_with_modified_sequence(self):
        data = {
            "Batch": ["A", "A", "B", "B"],
            "Modified sequence": ["seq1", "seq1", "seq2", "seq2"],
            "Proteins": ["P1", "P1", "P2", "P2"],
            "Gene names": ["G1", "G1", "G2", "G2"],
            "Reporter intensity corrected 1": [10, 20, 30, 40],
            "Identification metadata 1": ["imputed;", "imputed;", "", ""],
        }
        df = pd.DataFrame(data)
        result = prep.sum_peptide_intensities(df, run_lfq=True)
        expected_data = {
            "Batch": ["A", "B"],
            "Modified sequence": ["seq1", "seq2"],
            "Proteins": ["P1", "P2"],
            "Gene names": ["G1", "G2"],
            "Reporter intensity corrected 1": [30, 70],
            "Identification metadata 1": ["imputed;", ""],
        }
        expected_df = pd.DataFrame(expected_data)
        pd.testing.assert_frame_equal(result, expected_df)

    def test_aggregate_intensities_tmt_with_modified_sequence(self):
        data = {
            "Batch": ["A", "A", "B", "B"],
            "Modified sequence": ["seq1", "seq1", "seq2", "seq2"],
            "Proteins": ["P1", "P1", "P2", "P2"],
            "Gene names": ["G1", "G1", "G2", "G2"],
            "Reporter intensity corrected 1": [10, 20, 30, 40],
            "Reporter intensity corrected 2": [10, 20, 30, 40],
            "Reporter intensity corrected 3": [10, 20, 30, 40],
            "Reporter intensity corrected 4": [10, 20, 30, 40],
            "Reporter intensity corrected 5": [10, 20, 30, 40],
            "Reporter intensity corrected 6": [10, 20, 30, 40],
            "Reporter intensity corrected 7": [10, 20, 30, 40],
            "Reporter intensity corrected 8": [10, 20, 30, 40],
            "Reporter intensity corrected 9": [10, 20, 30, 40],
            "Reporter intensity corrected 10": [10, 20, 30, 40],
            "Reporter intensity corrected 11": [10, 20, 30, 40],
            "Identification metadata 1": ["imputed;", "imputed;", "", ""],
            "Identification metadata 2": ["imputed;", "imputed;", "", ""],
            "Identification metadata 3": ["imputed;", "imputed;", "", ""],
            "Identification metadata 4": ["imputed;", "imputed;", "", ""],
            "Identification metadata 5": ["imputed;", "imputed;", "", ""],
            "Identification metadata 6": ["imputed;", "imputed;", "", ""],
            "Identification metadata 7": ["imputed;", "imputed;", "", ""],
            "Identification metadata 8": ["imputed;", "imputed;", "", ""],
            "Identification metadata 9": ["imputed;", "imputed;", "", ""],
            "Identification metadata 10": ["imputed;", "imputed;", "", ""],
            "Identification metadata 11": ["imputed;", "imputed;", "", ""],
        }
        df = pd.DataFrame(data)
        result = prep.sum_peptide_intensities(df, run_lfq=False)
        expected_data = {
            "Batch": ["A", "B"],
            "Modified sequence": ["seq1", "seq2"],
            "Proteins": ["P1", "P2"],
            "Gene names": ["G1", "G2"],
            "Reporter intensity corrected 1": [30, 70],
            "Reporter intensity corrected 2": [30, 70],
            "Reporter intensity corrected 3": [30, 70],
            "Reporter intensity corrected 4": [30, 70],
            "Reporter intensity corrected 5": [30, 70],
            "Reporter intensity corrected 6": [30, 70],
            "Reporter intensity corrected 7": [30, 70],
            "Reporter intensity corrected 8": [30, 70],
            "Reporter intensity corrected 9": [30, 70],
            "Reporter intensity corrected 10": [30, 70],
            "Reporter intensity corrected 11": [30, 70],
            "Identification metadata 1": ["imputed;", ""],
            "Identification metadata 2": ["imputed;", ""],
            "Identification metadata 3": ["imputed;", ""],
            "Identification metadata 4": ["imputed;", ""],
            "Identification metadata 5": ["imputed;", ""],
            "Identification metadata 6": ["imputed;", ""],
            "Identification metadata 7": ["imputed;", ""],
            "Identification metadata 8": ["imputed;", ""],
            "Identification metadata 9": ["imputed;", ""],
            "Identification metadata 10": ["imputed;", ""],
            "Identification metadata 11": ["imputed;", ""],
        }
        expected_df = pd.DataFrame(expected_data)
        pd.testing.assert_frame_equal(result, expected_df)


class TestLoadAndNormalize:
    # Successfully loads and concatenates data from multiple batches
    def test_load_and_normalize_successful_load_and_concat_fp(self, mocker):
        # Mock dependencies
        mock_data_loader = mocker.Mock()
        mock_data_loader.load_data.return_value = [
            pd.DataFrame(
                {
                    "Modified sequence": [1, 2],
                    "Potential contaminant": ["", ""],
                    "Reverse": ["", ""],
                    "Batch": ["Batch1", "Batch1"],
                    "Reporter intensity corrected 1": [1, 2],
                }
            ),
            pd.DataFrame(
                {
                    "Modified sequence": [3, 4],
                    "Potential contaminant": ["", ""],
                    "Reverse": ["", ""],
                    "Batch": ["Batch2", "Batch2"],
                    "Reporter intensity corrected 1": [1, 2],
                }
            ),
        ]
        mock_data_loader.median_centering_within_batch.return_value = (
            pd.DataFrame({"A": [1, 2, 3, 4]}),
            pd.DataFrame(),
        )
        mock_data_loader.median_centering_ms1.return_value = (
            pd.DataFrame({"A": [1, 2, 3, 4]}),
            pd.DataFrame(),
        )
        mock_data_loader.impute_ms1_intensity.return_value = pd.DataFrame(
            {"A": [1, 2, 3, 4]}
        )
        mock_data_loader.redistribute_ms1_intensity.return_value = pd.DataFrame(
            {"A": [1, 2, 3, 4]}
        )

        mocker.patch("pandas.Series.to_csv")
        mocker.patch("pandas.DataFrame.to_csv")

        results_folder = "test_results"
        sample_annotation_df = pd.DataFrame(
            {
                "Sample name": ["Sample1", "Sample2"],
                "Batch Name": ["Batch1", "Batch2"],
                "TMT Channel": [1, 2],
                "is_reference": [False, True],
            }
        )
        data_type = "fp"

        # Call the function
        result_df = prep.load_and_normalize(
            mock_data_loader,
            results_folder,
            sample_annotation_df,
            data_type,
        )

        # Assertions
        assert not result_df.empty
        assert len(result_df) == 4

    def test_load_and_normalize_successful_load_and_concat_pp(self, mocker):
        # Mock dependencies
        mock_data_loader = mocker.Mock()
        mock_data_loader.load_data.return_value = [
            pd.DataFrame(
                {
                    "Modified sequence": [1, 2],
                    "Potential contaminant": ["", ""],
                    "Modifications": ["", ""],
                    "Reverse": ["", ""],
                    "Batch": ["Batch1", "Batch1"],
                    "Reporter intensity corrected 1": [1, 2],
                }
            ),
            pd.DataFrame(
                {
                    "Modified sequence": [3, 4],
                    "Potential contaminant": ["", ""],
                    "Modifications": ["", ""],
                    "Reverse": ["", ""],
                    "Batch": ["Batch2", "Batch2"],
                    "Reporter intensity corrected 1": [1, 2],
                }
            ),
        ]
        mock_data_loader.median_centering_within_batch.return_value = (
            pd.DataFrame({"A": [1, 2, 3, 4]}),
            pd.DataFrame(),
        )
        mock_data_loader.median_centering_ms1.return_value = (
            pd.DataFrame({"A": [1, 2, 3, 4]}),
            pd.DataFrame(),
        )
        mock_data_loader.impute_ms1_intensity.return_value = pd.DataFrame(
            {"A": [1, 2, 3, 4]}
        )
        mock_data_loader.redistribute_ms1_intensity.return_value = pd.DataFrame(
            {"A": [1, 2, 3, 4]}
        )

        mocker.patch("pandas.Series.to_csv")
        mocker.patch("pandas.DataFrame.to_csv")

        results_folder = "test_results"
        sample_annotation_df = pd.DataFrame(
            {
                "Sample name": ["Sample1", "Sample2"],
                "Batch Name": ["Batch1", "Batch2"],
                "TMT Channel": [1, 2],
                "is_reference": [False, True],
            }
        )
        data_type = "pp"

        # Call the function
        result_df = prep.load_and_normalize(
            mock_data_loader,
            results_folder,
            sample_annotation_df,
            data_type,
        )

        # Assertions
        assert not result_df.empty
        assert len(result_df) == 4


class TestImputation:
    """
    we need minimum 1 channel for imputation
    # check that minimum is chosen when minimum and max/100 when that's the case:   DONE
    # check imputation status added:                                                DONE
    # check that only patient channels are imputed:                                 DONE
    # check that there are no 0s:
    # perhaps just some rows for each test making it clear what is tested and what breaks !

    """

    # TODO: use fixture with readable dataframe
    def test_impute_data(self, df_tmt11, exp_df_tmt11_imputed_with_status):
        test_df = prep.impute_data(df_tmt11)
        exp_result = exp_df_tmt11_imputed_with_status
        # hacky way to test only imputation values
        test_df = test_df.loc[:, test_df.columns.str.startswith("Reporter")]
        exp_result = exp_result.loc[:, exp_result.columns.str.startswith("Reporter")]
        pd.testing.assert_frame_equal(test_df, exp_result, atol=0.01)

    def test_impute_data_imputation_status(
        self, df_tmt11, exp_df_tmt11_imputed_with_status
    ):
        test_df = prep.impute_data(df_tmt11)
        exp_result = exp_df_tmt11_imputed_with_status

        # hacky way to test only imputation values
        test_df = test_df.loc[
            :,
            (
                test_df.columns.str.startswith("Reporter")
                | test_df.columns.str.startswith("Identification")
            ),
        ]
        pd.testing.assert_frame_equal(test_df, exp_result, atol=0.01)


def test_impute_data_missing_ref():
    test_df = pd.DataFrame(
        [
            (
                19022.0,
                14528.0,
                32250.0,
                23406.0,
                621154.0,
                9854.9,
                141170.00,
                71932.0,
                90506.0,
                57724.0,
                57168.0,
                "x",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
            ),
            (
                np.nan,
                343490.0,
                np.nan,
                np.nan,
                310470.0,
                336740.0,
                550.37,
                np.nan,
                72777.0,
                57348.0,
                51079.0,
                "x;",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
            ),
            (
                50365.0,
                18493.0,
                69062.0,
                np.nan,
                16222.0,
                5868.5,
                3980.80,
                5669.2,
                np.nan,
                np.nan,
                np.nan,
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
            ),
            (
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
            ),
        ],
        columns=[f"Reporter intensity corrected {i}" for i in range(1, 12)]
        + [f"Identification metadata {i}" for i in range(1, 12)],
    )
    # print(test_df)
    test_df = prep.impute_data(test_df)
    expect_result = pd.DataFrame(
        [
            (
                19022.0,
                14528.0,
                32250.0,
                23406.0,
                621154.0,
                9854.9,
                141170.00,
                71932.0,
                90506.0,
                57724.0,
                57168.0,
                "x",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
            ),
            (
                550.37,
                343490.0,
                550.37,
                550.37,
                310470.0,
                336740.0,
                550.37,
                550.37,
                72777.0,
                57348.0,
                51079.0,
                "x;imputed;",
                "",
                "imputed;",
                "imputed;",
                "",
                "",
                "",
                "imputed;",
                "",
                "",
                "",
            ),
            (
                50365.0,
                18493.0,
                69062.0,
                690.62,
                16222.0,
                5868.5,
                3980.80,
                5669.2,
                np.nan,
                np.nan,
                np.nan,
                "",
                "",
                "",
                "imputed;",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
            ),
            (
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
            ),
        ],
        columns=[f"Reporter intensity corrected {i}" for i in range(1, 12)]
        + [f"Identification metadata {i}" for i in range(1, 12)],
    )
    pd.testing.assert_frame_equal(test_df, expect_result, atol=0.01)


class TestMedianCentering:
    # TODO: make simpler input dataframe where we can easily see the median
    def test_median_centering(self):
        test_df = pd.DataFrame(
            [
                (1, 1, 1, 1, 1, 1, 1, 1, 4, 3, 2),
                (1, 2, 1, 2, 1, 1, 1, 1, 4, 3, 2),
                (2, 2, 2, 2, 2, 2, 2, 2, 5, 6, 3),
                (2, 2, 2, 3, 6, 2, 2, 2, 5, 6, 3),
            ]
        )
        # average of the column medians is 2.23
        test_df, _ = prep.median_centering(test_df)
        expect_result = pd.DataFrame(
            [
                [1.48, 1.11, 1.48, 1.11, 1.48, 1.48, 1.48, 1.48, 1.98, 1.48, 1.78],
                [1.48, 2.23, 1.48, 2.23, 1.48, 1.48, 1.48, 1.48, 1.98, 1.48, 1.78],
                [2.97, 2.23, 2.97, 2.23, 2.97, 2.97, 2.97, 2.97, 2.47, 2.97, 2.67],
                [2.97, 2.23, 2.97, 3.34, 8.91, 2.97, 2.97, 2.97, 2.47, 2.97, 2.67],
            ]
        )
        pd.testing.assert_frame_equal(test_df, expect_result, atol=0.01)

    def test_median_centering_with_zeroes(self):
        test_df = pd.DataFrame(
            [
                (0, 1, 1, 1, 6, 1, 1, 1, 4, 3, 2),
                (2, 2, 2, 2, 2, 2, 2, 2, 5, 6, 3),
                (1, 2, 1, 3, 1, 1, 1, 1, 4, 3, 2),
                (2, 2, 2, 2, 2, 2, 2, 2, 5, 6, 3),
            ]
        )
        # average of the column medians is 2.318
        test_df, _ = prep.median_centering(test_df)
        expect_result = pd.DataFrame(
            [
                [np.nan, 1.16, 1.55, 1.16, 6.95, 1.55, 1.55, 1.55, 2.06, 1.55, 1.85],
                [2.32, 2.32, 3.09, 2.32, 2.32, 3.09, 3.09, 3.09, 2.58, 3.09, 2.78],
                [1.16, 2.32, 1.55, 3.48, 1.16, 1.55, 1.55, 1.55, 2.06, 1.55, 1.85],
                [2.32, 2.32, 3.09, 2.32, 2.32, 3.09, 3.09, 3.09, 2.58, 3.09, 2.78],
            ]
        )
        pd.testing.assert_frame_equal(test_df, expect_result, atol=0.01)

    def test_median_centering_too_many_zeroes_in_row(self):
        test_df = pd.DataFrame(
            [
                (0, 0, 0, 0, 6, 1, 1, 1, 4, 3, 2),
                (2, 2, 2, 2, 2, 2, 2, 2, 5, 6, 3),
                (1, 2, 1, 3, 1, 1, 1, 1, 4, 3, 2),
                (2, 2, 2, 2, 2, 2, 2, 2, 5, 6, 3),
            ]
        )
        # average of the column medians is 2.318
        test_df, _ = prep.median_centering(test_df)
        expect_result = pd.DataFrame(
            [
                [
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    8.18,
                    1.36,
                    1.36,
                    1.36,
                    2.18,
                    1.36,
                    1.82,
                ],
                [2.73, 2.73, 2.73, 2.73, 2.73, 2.73, 2.73, 2.73, 2.73, 2.73, 2.73],
                [1.36, 2.73, 1.36, 4.09, 1.36, 1.36, 1.36, 1.36, 2.18, 1.36, 1.82],
                [2.73, 2.73, 2.73, 2.73, 2.73, 2.73, 2.73, 2.73, 2.73, 2.73, 2.73],
            ]
        )
        pd.testing.assert_frame_equal(test_df, expect_result, atol=0.01)


#    def test_sample_correction(self):
#        # check that 0 turns into nan
#        test_df = pd.DataFrame(
#            [(1, 1, 1, 1, 6, 1, 1, 1, 4, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 5, 6, 3),
#             (1, 2, 1, 3, 1, 1, 1, 1, 4, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 5, 6, 3)])
#        test_df = preprocess.sample_correction(test_df)
#        print(test_df)
#        expect_result = pd.DataFrame(
#            [(1, 1, 1, 1, 6, 1, 1, 1, 4, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 5, 6, 3),
#             (1, 2, 1, 3, 1, 1, 1, 1, 4, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 5, 6, 3)])
#        pd.testing.assert_frame_equal(test_df, expect_result, atol=.01)


class TestRedistributeMS1:
    def test_ms1_imputation(self, long_df_tmt11):
        # TODO: how should we handle empty ref channels? sum or average?
        ref_channel_df = pd.DataFrame(
            {
                "TMT Channel": [10, 11, 10, 11, 10, 11, 10, 11],
                "Batch Name": [
                    "Batch_1",
                    "Batch_1",
                    "Batch_2",
                    "Batch_2",
                    "Batch_3",
                    "Batch_3",
                    "Batch_4",
                    "Batch_4",
                ],
            }
        )
        df = tmt_loader._impute_ms1_intensity(long_df_tmt11, ref_channel_df)
        # print(df)
        assert df.iloc[1]["MS1"] == 5.0  # missing MS1, same channel intensities
        assert (
            df.iloc[2]["MS1"] == 5.0
        )  # missing MS1, doubled all channel intensities ==> should not have any influence
        assert (
            df.iloc[5]["MS1"] == 6.0
        )  # missing MS1, average of other MS1s if QC channels have the same ratio relative to the summed
        # intensity
        assert (
            df.iloc[6]["MS1"] == 3.0
        )  # missing MS1, doubled intensity in reference channel, sum of all channels unchanged
        assert np.isnan(
            df.iloc[8]["MS1"]
        )  # missing MS1, no other batches to transfer from

    def test_redistribute_ms1_intensity(self, long_df_tmt11):
        ref_channel_df = pd.DataFrame(
            {
                "TMT Channel": [10, 11, 10, 11, 10, 11, 10, 11],
                "Batch Name": [
                    "Batch_1",
                    "Batch_1",
                    "Batch_2",
                    "Batch_2",
                    "Batch_3",
                    "Batch_3",
                    "Batch_4",
                    "Batch_4",
                ],
            }
        )
        df = tmt_loader._impute_ms1_intensity(long_df_tmt11, ref_channel_df)
        df = tmt_loader._redistribute_ms1_intensity(df)
        # print(df)
        assert (
            df.iloc[0]["Reporter intensity corrected 1"] == 0.25
        )  # 1.0 (TMT) / 20.0 (TMT sum) * 5.0 (MS1) = 0.25
        assert (
            df.iloc[1]["Reporter intensity corrected 1"] == 0.25
        )  # 1.0 / 20.0 * 5.0 = 0.25
        assert (
            df.iloc[2]["Reporter intensity corrected 1"] == 0.25
        )  # 1.0 / 20.0 * 5.0 = 0.25
        assert np.isnan(
            df.iloc[5]["Reporter intensity corrected 2"]
        )  # do not change abundance if TMT abundance is missing
        assert (
            df.iloc[8]["Reporter intensity corrected 1"] == 3.00
        )  # do not change abundance if MS1 is missing


class TestConvertLongToWideFormat:
    def test_convert_long_to_wide_format(self, long_df):
        # print(long_df)
        wide_df = prep.convert_long_to_wide_format(long_df)
        # print(wide_df)
        assert len(wide_df.index) == 4
        # TODO: check if this is the behavior we want with the duplicate Gene names and Proteins columns
        assert (
            wide_df.columns
            == [
                "Modified sequence",
                "Reporter intensity corrected 1 Batch_1",
                "Reporter intensity corrected 2 Batch_1",
                "Reporter intensity corrected 3 Batch_1",
                "Gene names Batch_1",
                "Proteins Batch_1",
                "Reporter intensity corrected 1 Batch_2",
                "Reporter intensity corrected 2 Batch_2",
                "Reporter intensity corrected 3 Batch_2",
                "Gene names Batch_2",
                "Proteins Batch_2",
            ]
        ).all()
        assert len(wide_df.dropna().index) == 2
        # TODO: add more assertions to check for gene name and protein aggregation

    def test_convert_long_to_wide_format_with_metadata(self, long_df):
        # print(long_df)
        long_df["Identification metadata 1"] = ""
        long_df["Identification metadata 2"] = ""
        long_df["Identification metadata 3"] = ""
        long_df["Irrelevant metadata 1"] = ""
        wide_df = prep.convert_long_to_wide_format(long_df, has_metadata_cols=True)
        # print(wide_df)
        assert (
            wide_df.columns
            == [
                "Modified sequence",
                "Reporter intensity corrected 1 Batch_1",
                "Reporter intensity corrected 2 Batch_1",
                "Reporter intensity corrected 3 Batch_1",
                "Identification metadata 1 Batch_1",
                "Identification metadata 2 Batch_1",
                "Identification metadata 3 Batch_1",
                "Reporter intensity corrected 1 Batch_2",
                "Reporter intensity corrected 2 Batch_2",
                "Reporter intensity corrected 3 Batch_2",
                "Identification metadata 1 Batch_2",
                "Identification metadata 2 Batch_2",
                "Identification metadata 3 Batch_2",
                "Gene names",
                "Proteins",
            ]
        ).all()
        assert len(wide_df.index) == 4
        assert len(wide_df.dropna().index) == 2


class TestDropDuplicateIndices:
    def test_removes_duplicate_indices(self):
        data = {"A": [1, 2, 3, 4]}
        index = [0, 1, 1, 2]
        df = pd.DataFrame(data, index=index)
        result = prep.drop_duplicate_indices(df)
        expected_data = {"A": [1, 2, 4]}
        expected_index = [0, 1, 2]
        expected_df = pd.DataFrame(expected_data, index=expected_index)
        pd.testing.assert_frame_equal(result, expected_df)


class TestMergeMs1Columns:
    def test_merge_ms1_columns_multiple_batches(self):
        data = {
            "Batch": ["Batch1", "Batch1", "Batch2", "Batch2"],
            "Modified sequence": ["Seq1", "Seq2", "Seq1", "Seq3"],
            "Intensity": [100, 200, 150, 300],
        }
        df = pd.DataFrame(data)
        result = prep.merge_ms1_columns(df)

        expected_data = {"Batch1": [100, 200, None], "Batch2": [150, None, 300]}
        expected_index = pd.Series(["Seq1", "Seq2", "Seq3"], name="Modified sequence")
        expected_df = pd.DataFrame(expected_data, index=expected_index)

        pd.testing.assert_frame_equal(result, expected_df)


class TestRenameColumnsWithSampleIds:
    def test_rename_columns_correctly(self):
        df = pd.DataFrame(
            {
                "Reporter intensity corrected 1 batch1": [10, 20],
                "Reporter intensity corrected 2 batch1": [30, 40],
                "Other column": [50, 60],
            }
        )
        channel_to_sample_id_dict = {
            "Reporter intensity corrected 1 batch1": "Sample1",
            "Reporter intensity corrected 2 batch1": "Sample2",
        }
        index_cols = ["Other column"]
        result_df = sample_mapping.rename_columns_with_sample_ids(
            df, channel_to_sample_id_dict, index_cols
        )
        expected_columns = ["Other column", "pat_Sample1", "pat_Sample2"]
        assert list(result_df.columns) == expected_columns


class TestGetDataLocation:
    # Returns list of file paths when valid directory and file type are provided
    def test_returns_list_of_file_paths(self, mocker):
        mocker.patch("os.path.isdir", return_value=True)
        mocker.patch(
            "os.walk",
            return_value=[
                ("/mocked_dir", ("subdir",), ("evidence.txt", "other.txt")),
                ("/mocked_dir/subdir", (), ("evidence.txt",)),
            ],
        )
        mocker.patch("topas_pipeline.preprocess_tools.is_valid_fp_file", return_value=True)
        mocker.patch("topas_pipeline.preprocess_tools.is_valid_pp_file", return_value=True)

        result = prep.get_data_location("/mocked_dir", "fp")
        assert result == ["/mocked_dir/evidence.txt", "/mocked_dir/subdir/evidence.txt"]


class TestIdentificationMetadata:
    """

    Check that there can already be imputed values in there before adding num_peptides: DONE
    Check that 0 does not give number of peptides: DONE
    Check that nan does not give number of peptides:
    Check that we get ; after:
    TODO: use fixture. add row with nan example

    Do we want to fix the xnum_peptides situations?
    """

    def test_mark_num_peptides(self):
        test_df = pd.DataFrame(
            [
                (
                    2,
                    0,
                    2,
                    2,
                    2,
                    2,
                    1,
                    2,
                    2,
                    2,
                    3,
                    "x",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                ),
                (
                    3,
                    1,
                    2,
                    2,
                    2,
                    2,
                    2,
                    2,
                    2,
                    2,
                    1,
                    "x;",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "y;",
                ),
            ],
            columns=[
                f"Unique peptides Reporter intensity corrected {i}"
                for i in range(1, 12)
            ]
            + [f"Identification metadata {i}" for i in range(1, 12)],
        )
        test_df = mark_num_peptides(test_df)
        print(test_df)
        expect_result = pd.DataFrame(
            [
                (
                    2,
                    0,
                    2,
                    2,
                    2,
                    2,
                    1,
                    2,
                    2,
                    2,
                    3,
                    "xnum_peptides=2;",
                    "",
                    "num_peptides=2;",
                    "num_peptides=2;",
                    "num_peptides=2;",
                    "num_peptides=2;",
                    "num_peptides=1;",
                    "num_peptides=2;",
                    "num_peptides=2;",
                    "num_peptides=2;",
                    "num_peptides=3;",
                ),
                (
                    3,
                    1,
                    2,
                    2,
                    2,
                    2,
                    2,
                    2,
                    2,
                    2,
                    1,
                    "x;num_peptides=3;",
                    "num_peptides=1;",
                    "num_peptides=2;",
                    "num_peptides=2;",
                    "num_peptides=2;",
                    "num_peptides=2;",
                    "num_peptides=2;",
                    "num_peptides=2;",
                    "num_peptides=2;",
                    "num_peptides=2;",
                    "y;num_peptides=1;",
                ),
            ],
            columns=[
                f"Unique peptides Reporter intensity corrected {i}"
                for i in range(1, 12)
            ]
            + [f"Identification metadata {i}" for i in range(1, 12)],
        )
        pd.testing.assert_frame_equal(test_df, expect_result, atol=0.01)

    def test_mark_num_peptides_with_na(self):
        test_df = pd.DataFrame(
            [
                (
                    2,
                    np.nan,
                    2,
                    2,
                    2,
                    2,
                    1,
                    2,
                    2,
                    2,
                    3,
                    "x",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                ),
                (
                    3,
                    1,
                    2,
                    2,
                    2,
                    2,
                    2,
                    2,
                    2,
                    2,
                    1,
                    "x;",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "y;",
                ),
            ],
            columns=[
                f"Unique peptides Reporter intensity corrected {i}"
                for i in range(1, 12)
            ]
            + [f"Identification metadata {i}" for i in range(1, 12)],
        )
        test_df = mark_num_peptides(test_df)
        print(test_df)
        expect_result = pd.DataFrame(
            [
                (
                    2,
                    np.nan,
                    2,
                    2,
                    2,
                    2,
                    1,
                    2,
                    2,
                    2,
                    3,
                    "xnum_peptides=2;",
                    "",
                    "num_peptides=2;",
                    "num_peptides=2;",
                    "num_peptides=2;",
                    "num_peptides=2;",
                    "num_peptides=1;",
                    "num_peptides=2;",
                    "num_peptides=2;",
                    "num_peptides=2;",
                    "num_peptides=3;",
                ),
                (
                    3,
                    1,
                    2,
                    2,
                    2,
                    2,
                    2,
                    2,
                    2,
                    2,
                    1,
                    "x;num_peptides=3;",
                    "num_peptides=1;",
                    "num_peptides=2;",
                    "num_peptides=2;",
                    "num_peptides=2;",
                    "num_peptides=2;",
                    "num_peptides=2;",
                    "num_peptides=2;",
                    "num_peptides=2;",
                    "num_peptides=2;",
                    "y;num_peptides=1;",
                ),
            ],
            columns=[
                f"Unique peptides Reporter intensity corrected {i}"
                for i in range(1, 12)
            ]
            + [f"Identification metadata {i}" for i in range(1, 12)],
        )
        pd.testing.assert_frame_equal(test_df, expect_result, atol=0.01)


class TestBatchHandling:
    """ """

    def test_get_unique_batches(self):
        batches = pd.DataFrame(
            [
                (1, "Sarcoma"),
                (1, "Sarcoma"),
                (1, "Sarcoma"),
                (6, "Sarcoma"),
                (1, "Chordoma"),
                (1, "Chordoma"),
                (1, "Chordoma"),
                (4, "Chordoma"),
            ]
        )
        batches.columns = ["Batch Name", "Cohort"]
        assert get_unique_batches(batches) == [
            [1, "Sarcoma"],
            [6, "Sarcoma"],
            [1, "Chordoma"],
            [4, "Chordoma"],
        ]

    def test_filter_evidence_files(self):
        evidence_files = [
            "/my/path/Sarcoma/Batch1_FP_Blabla/combined/txt/evidence.txt",
            "/my/path/Sarcoma/Batch2_FP_Blabla/combined/txt/evidence.txt",
            "/my/path/Sarcoma/Batch3_FP_Blabla/combined/txt/evidence.txt",
            "/my/path/Sarcoma/Batch4_FP_Blabla/combined/txt/evidence.txt",
            "/my/path/Sarcoma/Batch5_FP_Blabla/combined/txt/evidence.txt",
            "/my/path/Sarcoma/Batch6_FP_Blabla/combined/txt/evidence.txt",
            "/my/path/Chordoma/Batch1_FP_Blabla/combined/txt/evidence.txt",
            "/my/path/Chordoma/Batch4_FP_Blabla/combined/txt/evidence.txt",
        ]
        batches = [[1, "Sarcoma"], [6, "Sarcoma"], [1, "Chordoma"], [4, "Chordoma"]]
        filtered_evidence_files = prep.filter_evidence_files(
            evidence_files, "FP", batches
        )
        assert filtered_evidence_files == [
            "/my/path/Sarcoma/Batch1_FP_Blabla/combined/txt/evidence.txt",
            "/my/path/Sarcoma/Batch6_FP_Blabla/combined/txt/evidence.txt",
            "/my/path/Chordoma/Batch1_FP_Blabla/combined/txt/evidence.txt",
            "/my/path/Chordoma/Batch4_FP_Blabla/combined/txt/evidence.txt",
        ]

    def test_extract_batch_name(self):
        evidence_files = [
            "/my/path/Sarcoma/Batch1_FP_Blabla/combined/txt/evidence.txt",
            "/my/path/Sarcoma/Batch6_FP_Blabla/combined/txt/evidence.txt",
            "/my/path/Chordoma/Batch1_FP_Blabla/combined/txt/evidence.txt",
            "/my/path/Chordoma/Batch14_FP_Blabla/combined/txt/evidence.txt",
        ]

        extracted_batches = [data_loader.extract_batch_name(e) for e in evidence_files]
        assert extracted_batches == [
            "Sarcoma_Batch1",
            "Sarcoma_Batch6",
            "Chordoma_Batch1",
            "Chordoma_Batch14",
        ]

    def test_extract_experiment_name(self):
        evidence_files = [
            "/my/path/Sarcoma/Batch1_FP_Blabla/combined/txt/evidence.txt",
            "/my/path/Sarcoma/Batch6_FP_Blabla/combined/txt/evidence.txt",
            "/my/path/Chordoma/Batch1_FP_Blabla2/combined/txt/evidence.txt",
            "/my/path/Chordoma/Batch14_FP_Blabla/combined/txt/evidence.txt",
        ]

        extracted_batches = [
            data_loader.extract_experiment_name(e) for e in evidence_files
        ]
        assert extracted_batches == [
            "Batch1_FP_Blabla",
            "Batch6_FP_Blabla",
            "Batch1_FP_Blabla2",
            "Batch14_FP_Blabla",
        ]


# Creating dataframes from strings: https://towardsdatascience.com/67b0c2b71e6a
@pytest.fixture
def df_tmt11():
    df_string = """Sequence,                         Modified sequence,    Gene names,      Proteins, Reporter intensity corrected 1, Reporter intensity corrected 2, Reporter intensity corrected 3, Reporter intensity corrected 4, Reporter intensity corrected 5, Reporter intensity corrected 6, Reporter intensity corrected 7, Reporter intensity corrected 8, Reporter intensity corrected 9, Reporter intensity corrected 10, Reporter intensity corrected 11
DSDSWDADAFSVEDPVRK, DS(Phospho (STY))DS(Phospho (STY))WDADAFSVEDPVRK, GENE_A;GENE_B, PROT_A;PROT_B, , , , , , , , , , , 
DSDSWDADAFSVEDPVRK, DS(Phospho (STY))DS(Phospho (STY))WDADAFSVEDPVRK, GENE_B;GENE_A, PROT_B;PROT_A, 1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 1.0, 1.0
DSDSWDADAFSVEDPVRK, DS(Phospho (STY))DS(Phospho (STY))WDADAFSVEDPVRK, GENE_B;GENE_A, PROT_B;PROT_A, , , , , , , , , 6.0, 2.0, 2.0
     SSPTPESPTMLTK,      SS(Phospho (STY))PTPES(Phospho (STY))PTMLTK, GENE_C;GENE_B, PROT_A;PROT_B, 3.0, , , , , , , , , ,
     SSPTPESPTMLTK,      SS(Phospho (STY))PTPES(Phospho (STY))PTMLTK,        GENE_B,        PROT_B, 3.0, 3.0,    , 3.0, 3.0,    , 1.0, 2.0,    , 3.0, 2.0
     SSPTPESPTMLTK,      SS(Phospho (STY))PTPES(Phospho (STY))PTMLTK,        GENE_B,        PROT_B, 3.0,    , 3.0, 3.0,    , 3.0, 1.0, 2.0,    , 3.0, 2.0
     SSPTPESPTMLTK,      SS(Phospho (STY))PTPES(Phospho (STY))PTMLTK,        GENE_B,        PROT_B, 0.03,    , 3.0, 3.0,    , 3.0, 1.0, 1.0,    , 6.0,    
     QSPTPESPTMLTK,      QS(Phospho (STY))PTPES(Phospho (STY))PTMLTK, GENE_A;GENE_B,        PROT_B, , , , , , , , , 1.0, 3.0, 2.0
     KSPTPESPTMLTK,      KS(Phospho (STY))PTPES(Phospho (STY))PTMLTK, GENE_A;GENE_B,        PROT_B, 3.0, 2.0, , , , 1.0, 3.0, 2.0, , , 
"""
    df = pd.read_csv(io.StringIO(df_string), delimiter=",", skipinitialspace=True)
    for i in range(1, 12):
        df[f"Identification metadata {i}"] = ""
    return df


@pytest.fixture
def exp_df_tmt11_imputed_with_status():
    df_string = """ Reporter intensity corrected 1, Reporter intensity corrected 2, Reporter intensity corrected 3, Reporter intensity corrected 4, Reporter intensity corrected 5, Reporter intensity corrected 6, Reporter intensity corrected 7, Reporter intensity corrected 8, Reporter intensity corrected 9, Reporter intensity corrected 10, Reporter intensity corrected 11, Identification metadata 1, Identification metadata 2, Identification metadata 3, Identification metadata 4, Identification metadata 5, Identification metadata 6, Identification metadata 7, Identification metadata 8, Identification metadata 9, Identification metadata 10, Identification metadata 11
NaN,  NaN,   NaN,  NaN,  NaN,  NaN,  NaN,  NaN,  NaN,  NaN, NaN, '', '', '', '', '', '', '', '', '', '', ''
1.00, 2.00, 3.00, 1.00, 2.00, 3.00, 1.00, 2.00,  3.0,  1.0, 1.0, '', '', '', '', '', '', '', '', '', '', ''
0.06, 0.06, 0.06, 0.06, 0.06, 0.06, 0.06, 0.06,  6.0,  2.0, 2.0, 'imputed;', 'imputed;', 'imputed;', 'imputed;', 'imputed;', 'imputed;', 'imputed;', 'imputed;', '', '', ''
3.00, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03,  NaN,  NaN, NaN, '', 'imputed;', 'imputed;', 'imputed;', 'imputed;', 'imputed;', 'imputed;', 'imputed;', '', '', ''
3.00, 3.00, 0.03, 3.00, 3.00, 0.03, 1.00, 2.00,  NaN,  3.0, 2.0, '', '', 'imputed;', '', '', 'imputed;', '', '', '', '', ''
3.00, 0.03, 3.00, 3.00, 0.03, 3.00, 1.00, 2.00,  NaN,  3.0, 2.0, '', 'imputed;', '', '', 'imputed;', '', '', '', '', '', ''
0.03, 0.03, 3.00, 3.00, 0.03, 3.00, 1.00, 1.00,  NaN,  6.0, NaN, '', 'imputed;', '', '', 'imputed;', '', '', '', '', '', ''
0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03,  1.0,  3.0, 2.0, 'imputed;', 'imputed;', 'imputed;', 'imputed;', 'imputed;', 'imputed;', 'imputed;', 'imputed;', '', '', ''
3.00, 2.00, 0.03, 0.03, 0.03, 1.00, 3.00, 2.00,  NaN,  NaN, NaN, '', '', 'imputed;', 'imputed;', 'imputed;', '', '', '', '', '', ''
    """
    df = pd.read_csv(
        io.StringIO(df_string),
        dtype={f"Reporter intensity corrected {i}": float for i in range(1, 12)},
        delimiter=",",
        skipinitialspace=True,
        quotechar="'",
        keep_default_na=False,
        na_values="NaN",
    )
    return df


@pytest.fixture
def long_df():
    df_string = """Batch, Sequence,                         Modified sequence,    Gene names,      Proteins, Reporter intensity corrected 1, Reporter intensity corrected 2, Reporter intensity corrected 3
Batch_1, DSDSWDADAFSVEDPVRK, DS(Phospho (STY))DS(Phospho (STY))WDADAFSVEDPVRK, GENE_A;GENE_B, PROT_A;PROT_B, 0.5, 5.3, 9.2
Batch_1,      SSPTPESPTMLTK,      SS(Phospho (STY))PTPES(Phospho (STY))PTMLTK, GENE_C;GENE_B, PROT_A;PROT_B, 0.4, 7.3, 9.3
Batch_1, WSDSWDADAFSVEDPVRK, WS(Phospho (STY))DSWDADAFS(Phospho (STY))VEDPVRK, GENE_D;GENE_B, PROT_A;PROT_B, 0.3, 1.3, 9.4
Batch_2, DSDSWDADAFSVEDPVRK, DS(Phospho (STY))DS(Phospho (STY))WDADAFSVEDPVRK, GENE_B;GENE_A, PROT_B;PROT_A, 0.6, 2.3, 9.6
Batch_2,      SSPTPESPTMLTK,      SS(Phospho (STY))PTPES(Phospho (STY))PTMLTK,        GENE_B,        PROT_B, 0.8, 1.3, 9.8
Batch_2,      TSPTPESPTMLTK,      TS(Phospho (STY))PTPES(Phospho (STY))PTMLTK, GENE_A;GENE_B,        PROT_B, 0.7, 9.3, 9.7"""
    df = pd.read_csv(io.StringIO(df_string), delimiter=",", skipinitialspace=True)
    df["Modifications"] = ""
    df["Protein Names"] = ""
    df["Charge"] = 2
    df["m/z"] = 500.0
    df["Mass"] = 1000.0
    df["Missed cleavages"] = 0
    df["Length"] = 1
    df["Reverse"] = ""
    for i in range(1, 4):
        df[f"Identification metadata {i}"] = ""
    return df


@pytest.fixture
def long_df_tmt11():
    df_string = """Batch, Sequence,                         Modified sequence,    Gene names,      Proteins, MS1, Reporter intensity corrected 1, Reporter intensity corrected 2, Reporter intensity corrected 3, Reporter intensity corrected 4, Reporter intensity corrected 5, Reporter intensity corrected 6, Reporter intensity corrected 7, Reporter intensity corrected 8, Reporter intensity corrected 9, Reporter intensity corrected 10, Reporter intensity corrected 11
Batch_1, DSDSWDADAFSVEDPVRK, DS(Phospho (STY))DS(Phospho (STY))WDADAFSVEDPVRK, GENE_A;GENE_B, PROT_A;PROT_B, 5.0, 1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 1.0, 1.0
Batch_2, DSDSWDADAFSVEDPVRK, DS(Phospho (STY))DS(Phospho (STY))WDADAFSVEDPVRK, GENE_B;GENE_A, PROT_B;PROT_A,    , 1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 1.0, 1.0
Batch_2, DSDSWDADAFSVEDPVRK, DS(Phospho (STY))DS(Phospho (STY))WDADAFSVEDPVRK, GENE_B;GENE_A, PROT_B;PROT_A,    , 2.0, 4.0, 6.0, 2.0, 4.0, 6.0, 2.0, 4.0, 6.0, 2.0, 2.0
Batch_1,      SSPTPESPTMLTK,      SS(Phospho (STY))PTPES(Phospho (STY))PTMLTK, GENE_C;GENE_B, PROT_A;PROT_B, 4.0, 3.0, 2.0, 1.0, 3.0, 2.0, 1.0, 1.0, 1.0, 1.0, 4.0, 1.0   
Batch_2,      SSPTPESPTMLTK,      SS(Phospho (STY))PTPES(Phospho (STY))PTMLTK,        GENE_B,        PROT_B, 8.0, 3.0, 3.0,    , 3.0, 3.0,    , 1.0, 2.0,    , 3.0, 2.0
Batch_3,      SSPTPESPTMLTK,      SS(Phospho (STY))PTPES(Phospho (STY))PTMLTK,        GENE_B,        PROT_B,    , 3.0,    , 3.0, 3.0,    , 3.0, 1.0, 2.0,    , 3.0, 2.0
Batch_4,      SSPTPESPTMLTK,      SS(Phospho (STY))PTPES(Phospho (STY))PTMLTK,        GENE_B,        PROT_B,    , 3.0,    , 3.0, 3.0,    ,    ,    , 1.0,    , 6.0, 4.0   
Batch_2,      TSPTPESPTMLTK,      TS(Phospho (STY))PTPES(Phospho (STY))PTMLTK, GENE_A;GENE_B,        PROT_B, 0.7, 3.0, 2.0, 1.0, 3.0, 2.0, 1.0, 3.0, 2.0, 1.0, 3.0, 2.0
Batch_3,      QSPTPESPTMLTK,      QS(Phospho (STY))PTPES(Phospho (STY))PTMLTK, GENE_A;GENE_B,        PROT_B,    , 3.0, 2.0, 1.0, 3.0, 2.0, 1.0, 3.0, 2.0, 1.0, 3.0, 2.0
Batch_3,      KSPTPESPTMLTK,      KS(Phospho (STY))PTPES(Phospho (STY))PTMLTK, GENE_A;GENE_B,        PROT_B,    , 3.0, 2.0, 1.0, 3.0, 2.0, 1.0, 3.0, 2.0,    ,    , 2.0
"""
    df = pd.read_csv(io.StringIO(df_string), delimiter=",", skipinitialspace=True)
    df["Modifications"] = ""
    df["Protein Names"] = ""
    df["Charge"] = 2
    df["m/z"] = 500.0
    df["Mass"] = 1000.0
    df["Missed cleavages"] = 0
    df["Length"] = 1
    df["Reverse"] = ""
    for i in range(1, 4):
        df[f"Identification metadata {i}"] = ""
    return df
