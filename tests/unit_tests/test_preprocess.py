from pathlib import Path

import pytest
import pandas as pd

from topas_pipeline import preprocess
from topas_pipeline import config


class TestPreprocessRaw:
    # Successfully processes raw data with valid directories and annotations
    def test_successful_processing_with_valid_data(self, mocker):
        results_folder = "/path/to/results"
        sample_annotation_file = "/path/to/sample_annotation.csv"
        metadata_annotation = "/path/to/metadata_annotation.csv"
        run_simsi = True
        simsi_folder = "/path/to/simsi"
        data_types = ["type1", "type2"]

        mocker.patch("topas_pipeline.preprocess.preprocess_raw_data_type")
        mocker.patch("topas_pipeline.preprocess_tools.check_annot")

        preprocess.preprocess_raw(
            results_folder=results_folder,
            sample_annotation_file=sample_annotation_file,
            metadata_annotation=metadata_annotation,
            run_simsi=run_simsi,
            simsi_folder=simsi_folder,
            data_types=data_types,
        )

        assert preprocess.preprocess_raw_data_type.call_count == 2


class TestPreprocessRawDataType:
    def test_preprocessing_starts_when_no_preprocessed_file_exists(self, mocker):
        # Mocking os.path.exists to simulate no preprocessed file exists
        mocker.patch("os.path.exists", return_value=False)
        mock_logger = mocker.patch("topas_pipeline.preprocess.logger")

        # Creating a sample DataFrame for sample_annotation_df
        sample_annotation_df = pd.DataFrame(
            {"is_reference": [False, True], "Sample name": ["Sample1", "Sample2"]}
        )

        # Mocking other functions used within preprocess_raw_data_type
        mocker.patch(
            "topas_pipeline.preprocess.sample_annotation.filter_samples_by_metadata",
            return_value=sample_annotation_df,
        )
        mocker.patch(
            "topas_pipeline.preprocess.load_sample_data", return_value=pd.DataFrame()
        )
        mocker.patch(
            "topas_pipeline.preprocess.preprocess_fp", return_value=pd.DataFrame()
        )
        mocker.patch(
            "topas_pipeline.preprocess.sample_annotation.get_channel_to_sample_id_dict",
            return_value={},
        )
        mocker.patch("topas_pipeline.preprocess.utils.get_index_cols", return_value=[])
        mocker.patch(
            "topas_pipeline.preprocess.prep.rename_columns_with_sample_ids",
            return_value=pd.DataFrame(columns=["Gene names"]),
        )
        mocker.patch(
            "topas_pipeline.preprocess.get_prefix_renaming_dict", return_value={}
        )
        mocker.patch("pandas.DataFrame.to_csv")

        preprocess_config = config.Preprocessing(
            raw_data_location="raw_data",
            picked_fdr=0.01,
            fasta_file="fasta.fasta",
            fdr_num_threads=4,
            program="program",
            entity="entity",
            histologic_subtype="subtype",
            imputation=True,
            run_lfq=False,
            normalize_to_reference=False,
            debug=False,
        )
        preprocess.preprocess_raw_data_type(
            results_folder="results",
            run_simsi=True,
            simsi_folder="simsi",
            sample_annotation_df=sample_annotation_df,
            preprocessing_config=preprocess_config,
            data_type="fp",
        )

        # Assert that the logger was called with the start message
        mock_logger.info.assert_any_call("Preprocessing fp starts")


class TestGetPrefixRenamingDict:
    # Correctly renames columns with 'pat_' prefix when all patient_columns are present in df
    def test_correct_renaming(self):
        df = pd.DataFrame({"id": [1, 2], "name": ["Alice", "Bob"], "age": [25, 30]})
        patient_columns = ["id", "name"]
        expected_result = {"id": "pat_id", "name": "pat_name"}

        result = preprocess.get_prefix_renaming_dict(df, patient_columns)

        assert result == expected_result

        # Raises ValueError when a patient_column does not map correctly to its expected prefixed value

    def test_invalid_mapping_raises_value_error(self):
        df = pd.DataFrame({"id": [1, 2], "name": ["Alice", "Bob"], "age": [25, 30]})
        patient_columns = ["id", "name", "name", "age"]

        with pytest.raises(
            ValueError,
            match="Invalid mapping: name should map to pat_name, but got pat_age instead.",
        ):
            preprocess.get_prefix_renaming_dict(df, patient_columns)


class TestLoadSampleData:
    # Correctly loads and normalizes data when run_simsi is True
    def test_load_and_normalize_with_tmt(self, mocker):
        # Mock dependencies
        mock_get_evidence_files = mocker.patch(
            "topas_pipeline.preprocess_tools.get_evidence_files",
            return_value=["file1", "file2"],
        )
        mock_tmt_loader = mocker.patch("topas_pipeline.preprocess.TMTLoader")
        mock_simsi_loader = mocker.patch("topas_pipeline.preprocess.SimsiTMTLoader")
        mock_lfq_loader = mocker.patch("topas_pipeline.preprocess.LFQLoader")
        mock_load_and_normalize = mocker.patch(
            "topas_pipeline.preprocess_tools.load_and_normalize",
            return_value=pd.DataFrame(),
        )

        # Test data
        results_folder = Path("/path/to/results")
        sample_annotation_df = pd.DataFrame(
            {"sample": ["A", "B"], "is_reference": [False, True]}
        )
        simsi_folder = Path("/path/to/simsi")
        raw_data_location = Path("/path/to/raw")
        run_simsi = False
        run_lfq = False
        debug = False
        data_type = "type1"
        normalize_to_reference = False

        ref_channels = pd.DataFrame(
            {"sample": ["B"], "is_reference": [True]}, index=[1]
        )

        # Call the function
        df = preprocess.load_sample_data(
            results_folder,
            sample_annotation_df,
            run_simsi,
            simsi_folder,
            raw_data_location,
            run_lfq,
            debug,
            data_type,
            normalize_to_reference,
        )

        # Assertions
        mock_get_evidence_files.assert_called_once_with(
            sample_annotation_df, raw_data_location, data_type
        )
        mock_tmt_loader.assert_called_once_with(["file1", "file2"])
        assert not mock_simsi_loader.called
        assert not mock_lfq_loader.called
        mock_load_and_normalize.assert_called_once_with(
            mock_tmt_loader.return_value,
            results_folder,
            sample_annotation_df,
            data_type=data_type,
            normalize_to_reference=normalize_to_reference,
            debug=debug,
        )
        assert isinstance(df, pd.DataFrame)

    def test_load_and_normalize_with_simsi(self, mocker):
        # Mock dependencies
        mock_get_evidence_files = mocker.patch(
            "topas_pipeline.preprocess_tools.get_evidence_files",
            return_value=["file1", "file2"],
        )
        mock_simsi_loader = mocker.patch("topas_pipeline.preprocess.SimsiTMTLoader")
        mock_load_and_normalize = mocker.patch(
            "topas_pipeline.preprocess_tools.load_and_normalize",
            return_value=pd.DataFrame(),
        )

        # Test data
        results_folder = Path("/path/to/results")
        sample_annotation_df = pd.DataFrame(
            {"sample": ["A", "B"], "is_reference": [False, True]}
        )
        simsi_folder = Path("/path/to/simsi")
        raw_data_location = Path("/path/to/raw")
        run_simsi = True
        run_lfq = False
        debug = False
        data_type = "type1"
        normalize_to_reference = False

        # Call the function
        df = preprocess.load_sample_data(
            results_folder,
            sample_annotation_df,
            run_simsi,
            simsi_folder,
            raw_data_location,
            run_lfq,
            debug,
            data_type,
            normalize_to_reference,
        )

        # Assertions
        mock_get_evidence_files.assert_called_once_with(
            sample_annotation_df, raw_data_location, data_type
        )
        mock_simsi_loader.assert_called_once_with(
            ["file1", "file2"], results_folder, simsi_folder, data_type
        )
        mock_load_and_normalize.assert_called_once_with(
            mock_simsi_loader.return_value,
            results_folder,
            sample_annotation_df,
            data_type=data_type,
            normalize_to_reference=normalize_to_reference,
            debug=debug,
        )
        assert isinstance(df, pd.DataFrame)

    def test_load_and_normalize_with_lfq(self, mocker):
        # Mock dependencies
        mock_get_evidence_files = mocker.patch(
            "topas_pipeline.preprocess_tools.get_evidence_files",
            return_value=["file1", "file2"],
        )
        mock_lfq_loader = mocker.patch("topas_pipeline.preprocess.LFQLoader")
        mock_load_and_normalize = mocker.patch(
            "topas_pipeline.preprocess_tools.load_and_normalize",
            return_value=pd.DataFrame(),
        )

        # Test data
        results_folder = Path("/path/to/results")
        sample_annotation_df = pd.DataFrame(
            {"sample": ["A", "B"], "is_reference": [False, True]}
        )
        simsi_folder = Path("/path/to/simsi")
        raw_data_location = Path("/path/to/raw")
        run_simsi = False
        run_lfq = True
        debug = False
        data_type = "type1"
        normalize_to_reference = False

        # Call the function
        df = preprocess.load_sample_data(
            results_folder,
            sample_annotation_df,
            run_simsi,
            simsi_folder,
            raw_data_location,
            run_lfq,
            debug,
            data_type,
            normalize_to_reference,
        )

        # Assertions
        mock_get_evidence_files.assert_called_once_with(
            sample_annotation_df, raw_data_location, data_type
        )
        mock_lfq_loader.assert_called_once_with(["file1", "file2"])
        mock_load_and_normalize.assert_called_once_with(
            mock_lfq_loader.return_value,
            results_folder,
            sample_annotation_df,
            data_type=data_type,
            normalize_to_reference=normalize_to_reference,
            debug=debug,
        )
        assert isinstance(df, pd.DataFrame)


class TestPreprocessPp:
    # Test that the function runs through and all appropriate functions are called
    def test_function_runs_through_without_imputation(self, mocker):
        # Mocking the necessary functions with the correct import paths
        mock_remap_genes = mocker.patch(
            "topas_pipeline.picked_group.remap_gene_names",
            return_value=pd.DataFrame({"Gene names": ["gene1", "gene2"]}),
        )
        mock_create_metadata_columns = mocker.patch(
            "topas_pipeline.identification_metadata.create_metadata_columns",
            return_value=pd.DataFrame(),
        )
        mock_impute_data = mocker.patch(
            "topas_pipeline.preprocess_tools.impute_data", return_value=pd.DataFrame()
        )
        mock_filter_data = mocker.patch(
            "topas_pipeline.preprocess_tools.filter_data", return_value=pd.DataFrame()
        )
        mock_sum_peptide_intensities = mocker.patch(
            "topas_pipeline.preprocess_tools.sum_peptide_intensities",
            return_value=pd.DataFrame(),
        )
        mock_log_transform_intensities = mocker.patch(
            "topas_pipeline.preprocess_tools.log_transform_intensities",
            return_value=pd.DataFrame(),
        )
        mock_convert_long_to_wide = mocker.patch(
            "topas_pipeline.preprocess_tools.convert_long_to_wide_format",
            return_value=pd.DataFrame(),
        )
        mock_to_csv = mocker.patch("pandas.DataFrame.to_csv")

        df = pd.DataFrame({"uniprot_id": ["P12345", "P67890"]})
        fasta_file = "dummy.fasta"
        results_folder = "results"
        picked_fdr = 0.01
        fdr_num_threads = 4
        imputation = False
        debug = False
        run_lfq = False

        result_df = preprocess.preprocess_pp(
            df,
            results_folder,
            picked_fdr,
            fasta_file,
            fdr_num_threads,
            imputation,
            debug,
            run_lfq,
        )

        mock_remap_genes.assert_called_once()
        mock_create_metadata_columns.assert_called_once()
        assert not mock_impute_data.called
        assert not mock_to_csv.called
        mock_filter_data.assert_called_once()
        mock_sum_peptide_intensities.assert_called_once()
        mock_log_transform_intensities.assert_called_once()
        mock_convert_long_to_wide.assert_called_once()

    def test_function_runs_through_with_imputation(self, mocker):
        # Mocking the necessary functions with the correct import paths
        mock_remap_genes = mocker.patch(
            "topas_pipeline.picked_group.remap_gene_names",
            return_value=pd.DataFrame({"Gene names": ["gene1", "gene2"]}),
        )
        mock_create_metadata_columns = mocker.patch(
            "topas_pipeline.identification_metadata.create_metadata_columns",
            return_value=pd.DataFrame(),
        )
        mock_impute_data = mocker.patch(
            "topas_pipeline.preprocess_tools.impute_data", return_value=pd.DataFrame()
        )
        mock_filter_data = mocker.patch(
            "topas_pipeline.preprocess_tools.filter_data", return_value=pd.DataFrame()
        )
        mock_sum_peptide_intensities = mocker.patch(
            "topas_pipeline.preprocess_tools.sum_peptide_intensities",
            return_value=pd.DataFrame(),
        )
        mock_log_transform_intensities = mocker.patch(
            "topas_pipeline.preprocess_tools.log_transform_intensities",
            return_value=pd.DataFrame(),
        )
        mock_convert_long_to_wide = mocker.patch(
            "topas_pipeline.preprocess_tools.convert_long_to_wide_format",
            return_value=pd.DataFrame(),
        )
        mock_to_csv = mocker.patch("pandas.DataFrame.to_csv")

        df = pd.DataFrame({"uniprot_id": ["P12345", "P67890"]})
        fasta_file = "dummy.fasta"
        results_folder = "results"
        picked_fdr = 0.01
        fdr_num_threads = 4
        imputation = True
        debug = False
        run_lfq = False

        result_df = preprocess.preprocess_pp(
            df,
            results_folder,
            picked_fdr,
            fasta_file,
            fdr_num_threads,
            imputation,
            debug,
            run_lfq,
        )

        mock_remap_genes.assert_called_once()
        mock_create_metadata_columns.assert_called_once()
        mock_impute_data.assert_called_once()
        mock_to_csv.assert_called_once()
        mock_filter_data.assert_called_once()
        mock_sum_peptide_intensities.assert_called_once()
        mock_log_transform_intensities.assert_called_once()
        mock_convert_long_to_wide.assert_called_once()


class TestPreprocessFp:
    # Preprocesses DataFrame correctly with valid inputs
    def test_preprocess_fp_valid_inputs(self, mocker):
        import pandas as pd
        from topas_pipeline.preprocess import preprocess_fp

        # Mocking dependencies
        mock_picked_protein_grouping = mocker.patch(
            "topas_pipeline.picked_group.picked_protein_grouping",
            return_value=pd.DataFrame({"A": [1, 2], "B": [3, 4]}),
        )
        mock_filter_data = mocker.patch(
            "topas_pipeline.preprocess_tools.filter_data",
            return_value=pd.DataFrame({"A": [1, 2], "B": [3, 4]}),
        )
        mock_create_metadata_columns = mocker.patch(
            "topas_pipeline.identification_metadata.create_metadata_columns",
            return_value=pd.DataFrame({"A": [1, 2], "B": [3, 4]}),
        )
        mock_mark_detected_in_batch = mocker.patch(
            "topas_pipeline.identification_metadata.mark_detected_in_batch",
            return_value=pd.DataFrame({"A": [1, 2], "B": [3, 4]}),
        )
        mock_mark_num_peptides = mocker.patch(
            "topas_pipeline.identification_metadata.mark_num_peptides",
            return_value=pd.DataFrame({"A": [1, 2], "B": [3, 4]}),
        )
        mock_log_transform_intensities = mocker.patch(
            "topas_pipeline.preprocess_tools.log_transform_intensities",
            return_value=pd.DataFrame({"A": [1, 2], "B": [3, 4]}),
        )

        # Test data
        df = pd.DataFrame({"A": [1, 2], "B": [3, 4]})
        results_folder = "results"
        picked_fdr = 0.01
        fasta_file = "fasta_file.fasta"
        fdr_num_threads = 4
        imputation = True
        debug = False
        run_lfq = True

        # Call the function
        result_df = preprocess_fp(
            df,
            results_folder,
            picked_fdr,
            fasta_file,
            fdr_num_threads,
            imputation,
            debug,
            run_lfq,
        )

        # Assertions
        assert not result_df.empty
        assert "A" in result_df.columns
        assert "B" in result_df.columns

        mock_picked_protein_grouping.assert_called_once()
        mock_filter_data.assert_called_once()
        mock_create_metadata_columns.assert_called_once()
        mock_mark_detected_in_batch.assert_called_once()
        mock_mark_num_peptides.assert_called_once()
        mock_log_transform_intensities.assert_called_once()
