import pandas as pd
import pytest
from unittest import mock

from topas_pipeline.preprocess.file_processor.constants import (
    DIANN_REPORT_OUTPUT_COLUMNS,
    DIANN_PSM_OUTPUT_COLUMNS,
    MQ_OUTPUT_COLUMNS,
)
from topas_pipeline.preprocess.file_processor.diann_file_loader import (
    DiannFileLoader,
)


class TestDiannFileLoader:
    @pytest.fixture
    def diann_report_file_path(self, tmp_path):
        """Fixture to create a temporary DIANN report file for testing"""
        diann_report_content = (
            "Modified.Sequence\tPrecursor.Charge\tMs1.Normalised\tRun\tPEP\tRT\n"
            "DVFSGSDTDPDM(UniMod:35)AFC(UniMod:4)K\t2\t1234\tfile1_raw\t0.001\t5.00\n"
            "DVFSGSDTDPDMAFC(UniMod:4)K\t3\t1235\tfile2_raw\t0.002\t8.50\n"
        )
        tsv_path = tmp_path / "dia-quant-output" / "report.tsv"
        tsv_path.parent.mkdir(parents=True, exist_ok=True)
        tsv_path.write_text(diann_report_content)
        return str(tsv_path)

    @pytest.fixture
    def diann_psm_file_path(self, tmp_path):
        """Fixture to create a temporary DIANN PSM file for testing"""
        diann_psm_content = (
            "Modified Peptide\tPeptide\tCharge\tIs Decoy\tIs Contaminant\tHyperscore\tSpectrum\tAssigned Modifications\tProtein\tMapped Proteins\tGene\tMapped Genes\n"
            "DVFSGSDTDPDM[147]AFCK\tDVFSGSDTDPDMAFCK\t2\tfalse\tfalse\t123.45\txxx.100.100.0\t12M(15.9949), 15C(57.0214)\tsp|P12345|PROT_HUMAN\t\tGENE1\tGENE2, GENE3\n"
            "RGQTCVVHYTGM[147]LEDGKK\tRGQTCVVHYTGMLEDGKK\t3\tfalse\tfalse\t234.56\txxx.200.200.0\t12M(15.9949), 5C(57.0214)\tsp|P67890|PROT2_HUMAN\t\tGENE2\t\n"
            "RGQTCVVHYTGM[147]LEDGK\tRGQTCVVHYTGMLEDGK\t3\tfalse\ttrue\t234.56\txxx.300.300.0\t12M(15.9949), 5C(57.0214)\tsp|P67890|PROT3_HUMAN\t\tGENE3\t\n"
            "DVFSGSDTDPDMAFCK\tDVFSGSDTDPDMAFCK\t3\ttrue\tfalse\t234.56\txxx.400.400.0\t15C(57.0214)\trev_sp|P67890|PROT3_HUMAN\t\tGENE3\t\n"
        )
        tsv_path = tmp_path / "psm.tsv"
        tsv_path.write_text(diann_psm_content)
        return str(tsv_path)

    @mock.patch(
        "topas_pipeline.preprocess.file_processor.diann_file_loader.DiannFileLoader.extract_experiment_name_from_diann_report_path",
        return_value="TestExperiment",
    )
    def test_load_diann_report_file(
        self, mock_extract_experiment_name, diann_report_file_path
    ):
        """Test the load_diann_report_file method of DiannFileLoader"""
        loader = DiannFileLoader(diann_report_file_path)
        df = loader.load_diann_report_file(diann_report_file_path)

        # Check dataframe output structure and columns
        assert df.shape[0] == 2
        assert df.shape[1] == len(DIANN_REPORT_OUTPUT_COLUMNS)
        assert sorted(df.columns.tolist()) == sorted(DIANN_REPORT_OUTPUT_COLUMNS)

        # Check that the "Experiment" column is correctly populated
        pd.testing.assert_series_equal(
            df["Experiment"],
            pd.Series(["TestExperiment", "TestExperiment"], name="Experiment"),
        )

        # Check that the "Modified sequence" column is correctly processed
        pd.testing.assert_series_equal(
            df["Modified sequence"],
            pd.Series(
                ["_DVFSGSDTDPDM(Oxidation (M))AFCK_", "_DVFSGSDTDPDMAFCK_"],
                name="Modified sequence",
            ),
        )

        # Check that the Raw file column is correctly populated
        pd.testing.assert_series_equal(
            df["Raw file"], pd.Series(["file1_raw", "file2_raw"], name="Raw file")
        )

        # Check that the Type column is correctly populated
        pd.testing.assert_series_equal(
            df["Type"], pd.Series(["MULTI-MSMS", "MULTI-MSMS"], name="Type")
        )

        # Check that the Fraction column is correctly populated
        pd.testing.assert_series_equal(
            df["Fraction"], pd.Series([1, 1], name="Fraction")
        )

        # Check that the RT column is correctly populated
        pd.testing.assert_series_equal(
            df["RT"], pd.Series([300.00, 510.00], name="RT")
        )

    def test_load_diann_psm_file(self, diann_psm_file_path):
        """Test the load_diann_psm_file method of DiannFileLoader"""
        loader = DiannFileLoader("dummy_path")  # The path is not used in this test
        df = loader.load_psm_file(diann_psm_file_path)

        # Check dataframe output structure and columns
        assert df.shape[0] == 4
        assert df.shape[1] == len(DIANN_PSM_OUTPUT_COLUMNS)
        assert sorted(df.columns.tolist()) == sorted(DIANN_PSM_OUTPUT_COLUMNS)

        # Check that the "Modified sequence" column is correctly processed
        pd.testing.assert_series_equal(
            df["Modified sequence"],
            pd.Series(
                [
                    "_DVFSGSDTDPDM(Oxidation (M))AFCK_",
                    "_RGQTCVVHYTGM(Oxidation (M))LEDGKK_",
                    "_RGQTCVVHYTGM(Oxidation (M))LEDGK_",
                    "_DVFSGSDTDPDMAFCK_",
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

        # Check that the Leading proteins column is correctly populated
        pd.testing.assert_series_equal(
            df["Leading proteins"],
            pd.Series(
                [
                    "sp|P12345|PROT_HUMAN",
                    "sp|P67890|PROT2_HUMAN",
                    "sp|P67890|PROT3_HUMAN",
                    "rev_sp|P67890|PROT3_HUMAN",
                ],
                name="Leading proteins",
            ),
        )

        # Check that the Mapped Genes column is correctly populated
        pd.testing.assert_series_equal(
            df["Gene names"],
            pd.Series(
                ["GENE1;GENE2;GENE3", "GENE2", "GENE3", "GENE3"], name="Gene names"
            ),
        )

    @mock.patch(
        "topas_pipeline.preprocess.file_processor.diann_file_loader.DiannFileLoader.extract_experiment_name_from_diann_report_path",
        return_value="TestExperiment",
    )
    def test_load_single_file(
        self, mock_extract_experiment_name, diann_report_file_path, diann_psm_file_path
    ):
        """Test the load_single_file method of DiannFileLoader"""
        sample_annotation_df = pd.DataFrame(
            {
                "Sample name": ["Sample1"],
                "Cohort": ["TestCohort"],
            }
        )
        loader = DiannFileLoader(
            diann_report_file_path, sample_annotation_df=sample_annotation_df
        )
        df = loader.load_single_file(diann_report_file_path)

        # Check dataframe output structure and columns
        assert df.shape[0] == 2
        assert df.shape[1] == len(MQ_OUTPUT_COLUMNS) + 1 # TODO: RT to MQ_OUTPUT_COLUMNS?
        assert sorted(df.columns.tolist()) == sorted(MQ_OUTPUT_COLUMNS + ["RT"])

        # Check Modified sequence column of merged DataFrame
        pd.testing.assert_series_equal(
            df["Modified sequence"],
            pd.Series(
                ["_DVFSGSDTDPDM(Oxidation (M))AFCK_", "_DVFSGSDTDPDMAFCK_"],
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
            pd.Series([1234, 1235], name="Reporter intensity corrected 1"),
        )
