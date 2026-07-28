import pandas as pd
import pytest
from unittest import mock

from topas_pipeline.preprocess.file_processor.constants import (
    MQ_INPUT_COLUMNS_TMT,
    MQ_OUTPUT_COLUMNS_TMT,
)
from topas_pipeline.preprocess.file_processor.tmt_evidence_file_loader import (
    TMTEvidenceFileLoader,
)


class TestTMTEvidenceFileLoader:
    @pytest.fixture
    def evidence_file_path(self, tmp_path):
        """Fixture to create a temporary TMT evidence file for testing"""
        evidence_content = (
            "Modified sequence\tCharge\tProteins\tLeading proteins\tGene names\tIntensity\tRaw file\tid\tExperiment\tType\tMS/MS scan number\tPEP\tScore\tModifications\tReverse\tPotential contaminant\tFraction\tReporter intensity corrected 1\tReporter intensity corrected 2\tReporter intensity corrected 3\tReporter intensity corrected 4\tReporter intensity corrected 5\tReporter intensity corrected 6\tReporter intensity corrected 7\tReporter intensity corrected 8\tReporter intensity corrected 9\tReporter intensity corrected 10\tReporter intensity corrected 11\n"
            "_AAAAAAAATM(Oxidation (M))ALAAPSSPTPES(Phospho (STY))PTM(Oxidation (M))LTK_\t2\tsp|P12345|PROT_HUMAN\tsp|P12345|PROT_HUMAN\tGENE1\t1000\tfile1_raw\t1\t100\tMULTI=MSMS\t1234\t0.001\t33.211\tPhospho (STY)\t\t\t1\t1000\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n"
            "_AAAAAAAATM(Oxidation (M))ALAAPSS_\t3\tsp|P67890|PROT2_HUMAN\tsp|P67890|PROT2_HUMAN\tGENE2\t2000\tfile2_raw\t2\t100\tMULTI=MSMS\t1235\t0.002\t20.15\t\t\t1\t0\t2000\t0\t0\t0\t0\t0\t0\t0\t0\t0\n"
            "_BBBBBBBBBBK_\t2\tsp|P11111|PROT3_HUMAN\tsp|P11111|PROT3_HUMAN\tGENE3\t500\tfile3_raw\t3\t100\tMULTI=MSMS\t1236\t0.003\t15.0\t\t\t1\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n"
        )
        tsv_path = tmp_path / "evidence.txt"
        tsv_path.write_text(evidence_content)
        return str(tsv_path)

    @mock.patch(
        "topas_pipeline.preprocess.file_processor.tmt_evidence_file_loader.extract_batch_name",
        return_value="TestBatch",
    )
    def test_load_evidence_file(self, mock_extract_batch_name, evidence_file_path):
        """Test the load_evidence_file method of TMTEvidenceFileLoader"""
        loader = TMTEvidenceFileLoader(evidence_file_path)
        df = loader.load_evidence_file(evidence_file_path, usecols=MQ_INPUT_COLUMNS_TMT)

        # Check dataframe output structure and columns
        # The row with no measured reporter ions should be filtered out
        assert df.shape[0] == 2
        assert df.shape[1] == len(MQ_OUTPUT_COLUMNS_TMT)
        assert sorted(df.columns.tolist()) == sorted(MQ_OUTPUT_COLUMNS_TMT)

        # Check that the "Proteins" and "Leading proteins" columns are correctly processed
        pd.testing.assert_series_equal(
            df["Proteins"], pd.Series(["P12345", "P67890"], name="Proteins")
        )

        # Check that the Phosphosite notation is correctly converted in the "Modified sequence" column
        assert (
            df["Modified sequence"].iloc[0]
            == "_AAAAAAAATM(Oxidation (M))ALAAPSSPTPEpSPTM(Oxidation (M))LTK_"
        )

        # Check that columns requiring category dtype are correctly converted
        for column in ["Reverse", "Experiment", "Modifications", "Potential contaminant"]:
            assert df[column].dtype.name == "category"

        # Check that the "Batch" column is correctly populated based on the file path
        pd.testing.assert_series_equal(
            df["Batch"],
            pd.Series(["TestBatch", "TestBatch"], name="Batch"),
        )

    @mock.patch(
        "topas_pipeline.preprocess.file_processor.tmt_evidence_file_loader.extract_batch_name",
        return_value="TestBatch",
    )
    def test_load(self, mock_extract_batch_name, evidence_file_path):
        """Test the load method of TMTEvidenceFileLoader"""
        loader = TMTEvidenceFileLoader(evidence_file_path)
        all_batches = loader.load(usecols=MQ_INPUT_COLUMNS_TMT)

        # Check that the output is a list of dataframes
        assert isinstance(all_batches, list)
        assert len(all_batches) == 1
        assert isinstance(all_batches[0], pd.DataFrame)
