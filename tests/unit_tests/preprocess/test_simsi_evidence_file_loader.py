import pandas as pd
import pytest

from topas_pipeline.preprocess.file_processor.constants import (
    MQ_INPUT_COLUMNS_TMT,
    MQ_OUTPUT_COLUMNS_TMT,
)
from topas_pipeline.preprocess.file_processor.simsi_evidence_file_loader import (
    SimsiEvidenceFileLoader,
)

SIMSI_OUTPUT_COLUMNS = MQ_OUTPUT_COLUMNS_TMT + ["Transferred spectra count"]


class TestSimsiEvidenceFileLoader:
    @pytest.fixture
    def sample_annotation_df(self):
        """Fixture providing the batch/cohort mapping used to derive the "Batch" column"""
        return pd.DataFrame(
            {
                "Batch Name": ["1", "2"],
                "Cohort": ["Sarcoma", "Melanoma"],
            }
        )

    @pytest.fixture
    def evidence_file_path(self, tmp_path):
        """Fixture to create a temporary SIMSI evidence file for testing.

        SIMSI files use "Gene Names" (not "Gene names"), have no "Potential
        contaminant" column (it is derived from "Proteins"), and add a
        "Transferred spectra count" column.
        """
        columns = [
            "Modified sequence",
            "Charge",
            "Proteins",
            "Leading proteins",
            "Gene Names",
            "Intensity",
            "Raw file",
            "id",
            "Experiment",
            "Type",
            "MS/MS scan number",
            "PEP",
            "Score",
            "Modifications",
            "Reverse",
            "Fraction",
        ] + [f"Reporter intensity corrected {i}" for i in range(1, 12)] + [
            "Transferred spectra count"
        ]

        rows = [
            # phospho, batch1, kept
            [
                "_AAAAAAAATM(Oxidation (M))ALAAPSSPTPES(Phospho (STY))PTM(Oxidation (M))LTK_",
                2, "sp|P12345|PROT_HUMAN", "sp|P12345|PROT_HUMAN", "GENE1", 1000,
                "file1_raw", 1, "Batch1_Cohort", "MULTI=MSMS", 1234, 0.001, 33.211,
                "Phospho (STY)", "", 1,
            ] + [1000] + [0] * 10 + [5],
            # no phospho, batch1, kept
            [
                "_AAAAAAAATM(Oxidation (M))ALAAPSS_",
                3, "sp|P67890|PROT2_HUMAN", "sp|P67890|PROT2_HUMAN", "GENE2", 2000,
                "file2_raw", 2, "Batch1_Cohort", "MULTI=MSMS", 1235, 0.002, 20.15,
                "", "", 1,
            ] + [0, 2000] + [0] * 9 + [3],
            # contaminant, batch1, kept
            [
                "_CCCCCCCC_",
                2, "CON__Streptavidin", "CON__Streptavidin", "GENE3", 1500,
                "file3_raw", 3, "Batch1_Cohort", "MULTI=MSMS", 1236, 0.003, 10.5,
                "", "", 1,
            ] + [0, 0, 1500] + [0] * 8 + [1],
            # no reporter ions measured, batch1, filtered out
            [
                "_DDDDDDDD_",
                2, "sp|P99999|PROT4_HUMAN", "sp|P99999|PROT4_HUMAN", "GENE4", 300,
                "file4_raw", 4, "Batch1_Cohort", "MULTI=MSMS", 1237, 0.004, 5.0,
                "", "", 1,
            ] + [0] * 11 + [0],
            # batch2, kept
            [
                "_EEEEEEEE_",
                2, "sp|P55555|PROT5_HUMAN", "sp|P55555|PROT5_HUMAN", "GENE5", 800,
                "file5_raw", 5, "Batch2_Cohort", "MULTI=MSMS", 1238, 0.005, 12.0,
                "", "", 1,
            ] + [800] + [0] * 10 + [2],
        ]

        df = pd.DataFrame(rows, columns=columns)
        tsv_path = tmp_path / "simsi_evidence.txt"
        df.to_csv(tsv_path, sep="\t", index=False)
        return str(tsv_path)

    def test_load_evidence_file(
        self, evidence_file_path, sample_annotation_df
    ):
        """Test the load_evidence_file method of SimsiEvidenceFileLoader"""
        loader = SimsiEvidenceFileLoader(
            evidence_file_path, sample_annotation_df=sample_annotation_df
        )
        all_batches = loader.load_evidence_file(
            evidence_file_path, usecols=MQ_INPUT_COLUMNS_TMT
        )

        # load_evidence_file groups the output by "Batch", so it returns a
        # list of dataframes rather than a single dataframe
        assert isinstance(all_batches, list)
        assert len(all_batches) == 2
        df = pd.concat(all_batches, ignore_index=True)

        # Check dataframe output structure and columns
        # The row with no measured reporter ions should be filtered out
        assert df.shape[0] == 4
        assert df.shape[1] == len(SIMSI_OUTPUT_COLUMNS)
        assert sorted(df.columns.tolist()) == sorted(SIMSI_OUTPUT_COLUMNS)

        # Check that the "Proteins" and "Leading proteins" columns are correctly processed
        row1 = df.loc[df["Raw file"] == "file1_raw"].iloc[0]
        row2 = df.loc[df["Raw file"] == "file2_raw"].iloc[0]
        assert row1["Proteins"] == "P12345"
        assert row2["Proteins"] == "P67890"

        # Check that the Phosphosite notation is correctly converted in the "Modified sequence" column
        assert (
            row1["Modified sequence"]
            == "_AAAAAAAATM(Oxidation (M))ALAAPSSPTPEpSPTM(Oxidation (M))LTK_"
        )

        # Check that "Potential contaminant" is correctly derived from "Proteins"
        # containing "CON_", rather than being read directly from the file
        row3 = df.loc[df["Raw file"] == "file3_raw"].iloc[0]
        assert row3["Potential contaminant"] == "+"
        assert pd.isnull(row1["Potential contaminant"])

        # Check that the "Batch" column is correctly derived from the "Experiment"
        # column using the cohort/batch name mapping in sample_annotation_df
        assert row1["Batch"] == "Sarcoma_Batch1"
        row5 = df.loc[df["Raw file"] == "file5_raw"].iloc[0]
        assert row5["Batch"] == "Melanoma_Batch2"

        # Check that the row with no measured reporter ions was filtered out
        assert "file4_raw" not in df["Raw file"].values

    def test_load(self, evidence_file_path, sample_annotation_df):
        """Test the load method of SimsiEvidenceFileLoader"""
        loader = SimsiEvidenceFileLoader(
            evidence_file_path, sample_annotation_df=sample_annotation_df
        )
        all_batches = loader.load(usecols=MQ_INPUT_COLUMNS_TMT)

        # Check that the output is a list of dataframes, one per batch
        assert isinstance(all_batches, list)
        assert len(all_batches) == 2
        for batch_df in all_batches:
            assert isinstance(batch_df, pd.DataFrame)
            assert batch_df["Batch"].nunique() == 1
