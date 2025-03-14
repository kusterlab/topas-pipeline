import json

import pytest
import pandas as pd

import bin.clinical_process as cp
import bin.clinical_tools as clinical_tools


class TestClinicalProcess:
    def test_executes_clinical_process_data_type_for_each_data_type(self, mocker):
        def mock_clinical_process_data_type(*args, **kwargs):
            if "data_type" in kwargs:
                data_types_processed.append(kwargs["data_type"])

        data_types_processed = []
        data_types = ["fp", "pp"]

        mocker.patch(
            "bin.clinical_process.clinical_process_data_type",
            side_effect=mock_clinical_process_data_type,
        )
        cp.clinical_process(data_types=data_types)

        assert data_types_processed == data_types


class TestClinicalProcessDataType:
    # Processes data when all required files are present
    def test_processes_data_with_all_files_present(self, mocker):
        # Mocking dependencies
        mocker.patch("os.path.exists", return_value=False)
        mock_read_csv = mocker.patch("pandas.read_csv", return_value=pd.DataFrame())
        mocker.patch("bin.clinical_tools.phospho_annot", return_value=pd.DataFrame())
        mocker.patch("bin.clinical_tools.add_psp_urls", return_value=pd.DataFrame())
        mocker.patch(
            "bin.clinical_tools.prot_basket_annotation",
            return_value=(pd.DataFrame(), {}),
        )
        mocker.patch("json.dump")
        mocker.patch("pandas.DataFrame.to_csv")
        mocker.patch("builtins.open", mocker.mock_open())

        # Test data
        results_folder = "test_results"
        extra_kinase_annot = False
        debug = False
        prot_baskets = "test_baskets"
        pspFastaFile = "test_fasta"
        pspKinaseSubstrateFile = "test_kinase"
        pspAnnotationFile = "test_annotation"
        pspRegulatoryFile = "test_regulatory"
        data_type = "pp"

        # Call the function
        cp.clinical_process_data_type(
            results_folder,
            debug,
            prot_baskets,
            extra_kinase_annot,
            pspFastaFile,
            pspKinaseSubstrateFile,
            pspAnnotationFile,
            pspRegulatoryFile,
            data_type,
        )

        # # Assertions
        mock_read_csv.assert_called()
        clinical_tools.phospho_annot.assert_called()
        clinical_tools.add_psp_urls.assert_called()
        clinical_tools.prot_basket_annotation.assert_called()
        json.dump.assert_called()
        pd.DataFrame.to_csv.assert_called()


class TestMergeBasketsWithSubbaskets:
    # merge baskets and subbaskets correctly when both are non-empty
    def test_merge_non_empty_baskets_and_subbaskets(self):
        import pandas as pd
        from bin.clinical_process import merge_baskets_with_subbaskets

        data = {"basket": "fruit;vegetable", "sub_basket": "apple;carrot"}
        row = pd.Series(data)

        result = cp.merge_baskets_with_subbaskets(row)

        assert result == "fruit - apple;vegetable - carrot"


class TestGetUniqueBaskets:
    def test_handles_list_input_correctly(self):
        baskets = ['apple', 'banana', 'apple', 'orange']
        expected_output = 'apple;banana;orange'
        assert cp.get_unique_baskets(baskets) == expected_output
    
    def test_handles_string_input_correctly(self):
        baskets = 'apple;banana;apple;orange'
        expected_output = 'apple;banana;orange'
        assert cp.get_unique_baskets(baskets) == expected_output
    
        # Handles float input correctly
    def test_handles_float_input_correctly(self):
        baskets = 3.5
        expected_output = 3.5
        assert cp.get_unique_baskets(baskets) == expected_output


class TestReadAnnotationFiles:
        # Correctly reads annotation files from the specified results folder
    def test_reads_annotation_files_correctly(self, mocker):
        # Mocking the utils.get_index_cols function
        mocker.patch('bin.utils.get_index_cols', return_value=['index_col'])

        # Mocking the pd.read_csv function
        mock_annot = pd.DataFrame({'basket': ['basket1'], 'sub_basket': ['sub_basket1']})
        mocker.patch('pandas.read_csv', return_value=mock_annot)

        results_folder = 'test_folder'
        debug = False
        data_type = 'test_type'

        annot, annot_ref = cp.read_annotation_files(results_folder, debug, data_type)

        assert annot is not None
        assert annot_ref is None
        assert 'sub_basket' in annot.columns
        assert 'basket' in annot.columns