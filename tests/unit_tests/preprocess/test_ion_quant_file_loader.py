import pandas as pd
import pytest
from unittest import mock

from topas_pipeline.preprocess.file_processor.constants import (
    IONQUANT_COMBINED_ION_INPUT_COLUMNS,
    IONQUANT_COMBINED_ION_OUTPUT_COLUMNS,
    IONQUANT_PSM_OUTPUT_COLUMNS,
    MQ_OUTPUT_COLUMNS,
)
from topas_pipeline.preprocess.file_processor.ion_quant_file_loader import (
    IonQuantFileLoader,
)


class TestIonQuantFileLoader:
    @pytest.fixture
    def combined_ion_file_path(self, tmp_path):
        """Fixture to create a temporary IonQuant combined_ion file for testing"""
        combined_ion_content = (
            "Modified Sequence\tCharge\tProtein\tMapped Proteins\tGene\tMapped Genes\tfile1_raw Intensity\tfile2_raw Intensity\tfile1_raw Apex Retention Time\tfile2_raw Apex Retention Time\n"
            "AAAAAAAAAAATSGSGGC[57.0215]PPAPGLES[79.9663]GVGAVGC[57.0215]GY[79.9663]PR\t3\tsp|Q5VZB9|DMRTA_HUMAN\t\tDMRTA1\t\t0.0\t1659023.0\t305.0\t306.0\n"
            "AAAAAAAATMALAAPSSPTPESPT[79.9663]M[15.9949]LTK\t3\tsp|Q9NQS7|INCE_HUMAN\t\tINCENP\t\t0.0\t1480978.0\t1126.0\t1125.0\n"
        )
        tsv_path = tmp_path / "ionquant-output" / "combined_ion.tsv"
        tsv_path.parent.mkdir(parents=True, exist_ok=True)
        tsv_path.write_text(combined_ion_content)
        return str(tsv_path)

    @pytest.fixture
    def psm_file_path(self, tmp_path):
        """Fixture to create a temporary IonQuant PSM file for testing"""
        ionquant_psm_content = (
            "Modified Peptide\tPeptide\tCharge\tIs Decoy\tIs Contaminant\tHyperscore\tSpectrum\tAssigned Modifications\tProtein\tMapped Proteins\tGene\tMapped Genes\tProbability\n"
            "AAAAAAAAAAATSGSGGCPPAPGLES[167]GVGAVGCGY[243]PR\tAAAAAAAAAAATSGSGGCPPAPGLESGVGAVGCGYPR\t3\tfalse\tfalse\t123.45\txxx.100.100.0\t26S(79.9663), 35Y(79.9663)\tsp|Q5VZB9|DMRTA_HUMAN\t\tDMRTA1\t\t0.999\n"
            "RGQTCVVHYTGM[147]LEDGKK\tRGQTCVVHYTGMLEDGKK\t3\tfalse\tfalse\t234.56\txxx.200.200.0\t12M(15.9949), 5C(57.0214)\tsp|P67890|PROT2_HUMAN\t\tGENE2\t\t0.998\n"
            "RGQTCVVHYTGM[147]LEDGK\tRGQTCVVHYTGMLEDGK\t3\tfalse\ttrue\t234.56\txxx.300.300.0\t12M(15.9949), 5C(57.0214)\tsp|P67890|PROT3_HUMAN\t\tGENE3\t\t0.997\n"
            "AAAAAAAATMALAAPSSPTPESPT[181]M[147]LTK\tAAAAAAAATMALAAPSSPTPESPTMLTK\t3\ttrue\tfalse\t234.56\txxx.400.400.0\t24T(79.9663),25M(15.9949)\trev_sp|Q9NQS7|INCE_HUMAN\t\tINCENP\t\t0.5\n"
        )
        tsv_path = tmp_path / "ionquant-output" / "R1" / "psm.tsv"
        tsv_path.parent.mkdir(parents=True, exist_ok=True)
        tsv_path.write_text(ionquant_psm_content)
        return str(tsv_path)

    @mock.patch(
        "topas_pipeline.preprocess.file_processor.ion_quant_file_loader.IonQuantFileLoader.extract_experiment_name_from_ionquant_path",
        return_value="TestExperiment",
    )
    def test_load_ionquant_combined_ion_file(
        self, mock_extract_experiment_name, combined_ion_file_path
    ):
        """Test the load_ionquant_combined_ion_file method of IonQuantFileLoader"""
        loader = IonQuantFileLoader(combined_ion_file_path)
        df = loader.load_ionquant_combined_ion_file(combined_ion_file_path)

        # Check dataframe output structure and columns
        assert df.shape[0] == 2
        assert df.shape[1] == len(IONQUANT_COMBINED_ION_OUTPUT_COLUMNS)
        assert sorted(df.columns.tolist()) == sorted(
            IONQUANT_COMBINED_ION_OUTPUT_COLUMNS
        )

        # Check that the "Experiment" column is correctly populated
        pd.testing.assert_series_equal(
            df["Experiment"],
            pd.Series(["TestExperiment", "TestExperiment"], name="Experiment"),
        )

        # Check that the "Modified sequence" column is correctly processed
        pd.testing.assert_series_equal(
            df["Modified sequence"],
            pd.Series(
                [
                    "_AAAAAAAAAAATSGSGGCPPAPGLES(Phospho (STY))GVGAVGCGY(Phospho (STY))PR_",
                    "_AAAAAAAATMALAAPSSPTPESPT(Phospho (STY))M(Oxidation (M))LTK_",
                ],
                name="Modified sequence",
            ),
        )

        # Check that the Raw file column is correctly populated
        pd.testing.assert_series_equal(
            df["Raw file"],
            pd.Series(
                [
                    "file1_raw;file2_raw",
                    "file1_raw;file2_raw",
                ],
                name="Raw file",
            ),
        )

        # Check that the Type column is correctly populated
        pd.testing.assert_series_equal(
            df["Type"], pd.Series(["MULTI-MSMS", "MULTI-MSMS"], name="Type")
        )

        # Check that the Charge column is correctly populated
        pd.testing.assert_series_equal(df["Charge"], pd.Series([3, 3], name="Charge"))

        # Check that the RT column is correctly populated
        pd.testing.assert_series_equal(
            df["RT"], pd.Series([305.5, 1125.5], name="RT")
        )

    def test_load_ionquant_psm_file(self, psm_file_path):
        """Test the load_psm_files method of IonQuantFileLoader"""
        loader = IonQuantFileLoader("dummy_path")  # The path is not used in this test
        df = loader.load_psm_files([psm_file_path]).reset_index()

        # Check dataframe output structure and columns
        assert df.shape[0] == 4
        assert df.shape[1] == len(IONQUANT_PSM_OUTPUT_COLUMNS)
        assert sorted(df.columns.tolist()) == sorted(IONQUANT_PSM_OUTPUT_COLUMNS)

        # Check that the "Modified sequence" column is correctly processed
        pd.testing.assert_series_equal(
            df["Modified sequence"],
            pd.Series(
                [
                    "_AAAAAAAAAAATSGSGGCPPAPGLES(Phospho (STY))GVGAVGCGY(Phospho (STY))PR_",
                    "_RGQTCVVHYTGM(Oxidation (M))LEDGKK_",
                    "_RGQTCVVHYTGM(Oxidation (M))LEDGK_",
                    "_AAAAAAAATMALAAPSSPTPESPT(Phospho (STY))M(Oxidation (M))LTK_",
                ],
                name="Modified sequence",
            ),
        )

        # Check that the "Reverse" column is correctly populated
        pd.testing.assert_series_equal(
            df["Reverse"], pd.Series(["", "", "", "+"], name="Reverse")
        )

        # Check that the "Potential contaminant" column is correctly populated
        pd.testing.assert_series_equal(
            df["Potential contaminant"],
            pd.Series(["", "", "+", ""], name="Potential contaminant"),
        )

        # Check that the MS/MS scan number column is correctly populated
        pd.testing.assert_series_equal(
            df["MS/MS scan number"],
            pd.Series([100, 200, 300, 400], name="MS/MS scan number"),
        )

    @mock.patch(
        "topas_pipeline.preprocess.file_processor.ion_quant_file_loader.IonQuantFileLoader.extract_experiment_name_from_ionquant_path",
        return_value="TestExperiment",
    )
    def test_load_single_file(
        self, mock_extract_experiment_name, combined_ion_file_path, psm_file_path
    ):
        """Test the load_single_file method of IonQuantFileLoader"""
        sample_annotation_df = pd.DataFrame(
            {
                "Sample name": ["Sample1"],
                "Cohort": ["TestCohort"],
            }
        )
        loader = IonQuantFileLoader(
            combined_ion_file_path, sample_annotation_df=sample_annotation_df
        )
        df = loader.load_single_file(combined_ion_file_path)

        # Check dataframe output structure
        output_columns = [x for x in MQ_OUTPUT_COLUMNS if x != "Fraction"] + [
            "Leading gene",
            "RT"
        ]
        assert df.shape[0] == 2
        assert df.shape[1] == len(output_columns)
        assert sorted(df.columns.tolist()) == sorted(output_columns)

        # Check Modified sequence column of merged DataFrame
        pd.testing.assert_series_equal(
            df["Modified sequence"],
            pd.Series(
                [
                    "_AAAAAAAAAAATSGSGGCPPAPGLEpSGVGAVGCGpYPR_",
                    "_AAAAAAAATMALAAPSSPTPESPpTM(Oxidation (M))LTK_",
                ],
                name="Modified sequence",
            ),
        )

        # Check that the Batch column is correctly populated
        pd.testing.assert_series_equal(
            df["Batch"],
            pd.Series(
                ["TestCohort_BatchTestExperiment", "TestCohort_BatchTestExperiment"],
                name="Batch",
            ),
        )

        # Check that the Reporter intensity corrected 1 column is correctly populated
        pd.testing.assert_series_equal(
            df["Reporter intensity corrected 1"],
            pd.Series([829511.5, 740489.0], name="Reporter intensity corrected 1"),
        )
