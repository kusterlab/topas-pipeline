from pandas import Index, DataFrame

from topas_pipeline import sample_annotation as sa


class TestGetSampleQcLotMappingDf:
    """Test suite for get_sample_qc_lot_mapping_df function."""

    def test_basic_functionality(self):
        """Test basic functionality with normal inputs."""
        sample_columns = Index(
            [
                "Reporter intensity corrected 10 Cohort1_BatchA",
                "Reporter intensity corrected 11 Cohort1_BatchA",
                "Reporter intensity corrected 10 Cohort2_BatchB",
            ]
        )
        sample_annotation_df = DataFrame(
            {
                "Batch Name": ["A", "A", "B"],
                "TMT Channel": [10, 11, 10],
                "QC Lot": [1.0, 2.0, 1.0],
                "is_reference": [True, False, True],
            }
        )

        result = sa.get_sample_qc_lot_mapping_df(sample_columns, sample_annotation_df)

        assert isinstance(result, DataFrame)
        assert list(result.columns) == [
            "batch",
            "channel",
            "QC Lot",
            "is_reference",
            "QC Lot group",
        ]
        assert len(result) == 3
        assert (
            result.loc["Reporter intensity corrected 10 Cohort1_BatchA", "batch"] == "A"
        )
        assert (
            result.loc["Reporter intensity corrected 10 Cohort1_BatchA", "channel"]
            == 10
        )
        assert (
            result.loc["Reporter intensity corrected 10 Cohort1_BatchA", "QC Lot"]
            == 1.0
        )
