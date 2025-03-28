import pandas as pd

import topas_pipeline.clinical_annotation as cp
from topas_pipeline import clinical_tools
from topas_pipeline import config


class TestClinicalProcess:
    def test_executes_add_clinical_annotations_data_type_for_each_data_type(
        self, mocker
    ):
        def mock_add_clinical_annotations_data_type(*args, **kwargs):
            if "data_type" in kwargs:
                data_types_processed.append(kwargs["data_type"])

        data_types_processed = []
        data_types = ["fp", "pp"]

        mocker.patch(
            "topas_pipeline.clinical_annotation.add_clinical_annotations_data_type",
            side_effect=mock_add_clinical_annotations_data_type,
        )
        cp.add_clinical_annotations(data_types=data_types, debug=False)

        assert data_types_processed == data_types


class TestClinicalProcessDataType:
    # Processes data when all required files are present
    def test_processes_data_with_all_files_present(self, mocker):
        # Mocking dependencies
        mocker.patch("os.path.exists", return_value=False)
        mock_read_csv = mocker.patch("pandas.read_csv", return_value=pd.DataFrame())
        mocker.patch(
            "topas_pipeline.clinical_tools.add_phospho_annotations", return_value=pd.DataFrame()
        )
        mocker.patch(
            "topas_pipeline.clinical_tools.add_psp_urls", return_value=pd.DataFrame()
        )
        mocker.patch(
            "topas_pipeline.clinical_tools.add_topas_annotations",
            return_value=pd.DataFrame(),
        )
        mocker.patch("json.dump")
        mocker.patch("pandas.DataFrame.to_csv")
        mocker.patch("builtins.open", mocker.mock_open())

        # Test data
        results_folder = "test_results"
        clinic_proc_config = config.ClinicProc(
            prot_baskets="test_topas_annotation_file",
            pspFastaFile="test_fasta",
            pspKinaseSubstrateFile="test_kinase",
            pspAnnotationFile="test_annotation",
            pspRegulatoryFile="test_regulatory",
            extra_kinase_annot="",
        )
        data_type = "pp"

        # Call the function
        cp.add_clinical_annotations_data_type(
            results_folder,
            clinic_proc_config,
            data_type,
        )

        # # Assertions
        mock_read_csv.assert_called()
        clinical_tools.add_phospho_annotations.assert_called()
        clinical_tools.add_psp_urls.assert_called()
        clinical_tools.add_topas_annotations.assert_called()
        pd.DataFrame.to_csv.assert_called()


class TestReadAnnotationFiles:
    # Correctly reads annotation files from the specified results folder
    def test_reads_annotation_files_correctly(self, mocker):
        # Mocking the utils.get_index_cols function
        mocker.patch("topas_pipeline.utils.get_index_cols", return_value=["index_col"])

        # Mocking the pd.read_csv function
        mock_annot = pd.DataFrame(
            {
                "TOPAS_score": ["topas1"],
                "TOPAS_subscore": ["subtopas1"],
                "POI_category": ["Adaptor protein"],
            }
        )
        mocker.patch("pandas.read_csv", return_value=mock_annot)

        results_folder = "test_folder"
        data_type = "test_type"

        annot = cp.read_annotated_expression_file(results_folder, data_type)

        assert annot is not None
        assert "TOPAS_score" in annot.columns
        assert "TOPAS_subscore" in annot.columns
