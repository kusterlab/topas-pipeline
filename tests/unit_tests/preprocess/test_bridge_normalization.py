import pandas as pd
import numpy as np
import pytest

from topas_pipeline.preprocess import bridge_normalization as bn


class TestRowWiseNormalize:

    @pytest.fixture
    def mapping_single_lot(self) -> pd.DataFrame:
        return pd.DataFrame(
            {
                "batch": ["B1", "B1", "B2", "B2"],
                "QC Lot": ["L1", "L1", "L1", "L1"],
                "is_reference": [True, False, True, False],
            },
            index=["S1", "S2", "S3", "S4"],
        )

    def test_no_correction_when_refs_equal(self, mapping_single_lot: pd.DataFrame):
        # Both batches have identical reference intensity (log10 scale)
        row = pd.Series(
            [2.0, 1.0, 2.0, 1.0],
            index=mapping_single_lot.index,
            name="intensity corrected",
        )

        result = bn.row_wise_normalize(row, mapping_single_lot)

        pd.testing.assert_series_equal(result, row, check_names=False)

    def test_batch_shift_is_corrected(self, mapping_single_lot: pd.DataFrame):
        """
        B1 ref = 2.0  (real = 100)
        B2 ref = 1.0  (real = 10)

        mix median in real space:
        mean_realspace(B1)=2.0
        mean_realspace(B2)=1.0
        median = 1.5

        corrections:
        B1 = 1.5 - 2.0 = -0.5
        B2 = 1.5 - 1.0 = +0.5
        """
        row = pd.Series(
            [2.0, 0.5, 1.0, 0.2],
            index=mapping_single_lot.index,
        )

        expected = pd.Series(
            [1.5, 0.0, 1.5, 0.7],  # corrected values
            index=row.index,
            name="intensity corrected",
        )

        result = bn.row_wise_normalize(row, mapping_single_lot)

        pd.testing.assert_series_equal(result.round(6), expected.round(6))

    def test_multiple_qc_lots_same_batch(self):
        mapping = pd.DataFrame(
            {
                "batch": ["B1", "B1", "B1", "B1"],
                "QC Lot": ["L1", "L1", "L2", "L2"],
                "is_reference": [True, False, True, False],
            },
            index=["S1", "S2", "S3", "S4"],
        )

        """
        L1 ref = 2.0
        L2 ref = 1.0
        median = 1.5

        corrections:
        L1 = -0.5
        L2 = +0.5
        batch mean correction = 0
        """
        row = pd.Series(
            [2.0, 0.3, 1.0, 0.1], index=mapping.index, name="intensity corrected"
        )

        result = bn.row_wise_normalize(row, mapping)

        pd.testing.assert_series_equal(result.round(6), row.round(6))

    def test_nan_is_reference_defaults_false(self):
        mapping = pd.DataFrame(
            {
                "batch": ["B1", "B1", "B1", "B1"],
                "QC Lot": ["L1", "L1", "L2", "L2"],
                "is_reference": [True, np.nan, True, np.nan],
            },
            index=["S1", "S2", "S3", "S4"],
        )

        """
        L1 ref = 2.0
        L2 ref = 1.0
        median = 1.5

        corrections:
        L1 = -0.5
        L2 = +0.5
        batch mean correction = 0
        """
        row = pd.Series(
            [2.0, 0.3, 1.0, 0.1], index=mapping.index, name="intensity corrected"
        )

        result = bn.row_wise_normalize(row, mapping)

        pd.testing.assert_series_equal(result.round(6), row.round(6))

    def test_multiple_batches_multiple_qc_lots(self):
        """
        Setup:

        Batch B1:
            L1 ref = 2.0
            L2 ref = 1.0

        Batch B2:
            L1 ref = 1.0
            L2 ref = 0.0

        Per-lot medians:
            L1 median = median(2.0, 1.0) = 1.5
            L2 median = median(1.0, 0.0) = 0.5

        Corrections:
            B1-L1 = 1.5 - 2.0 = -0.5
            B1-L2 = 0.5 - 1.0 = -0.5
            → B1 batch correction = mean(-0.5, -0.5) = -0.5

            B2-L1 = 1.5 - 1.0 = +0.5
            B2-L2 = 0.5 - 0.0 = +0.5
            → B2 batch correction = mean(+0.5, +0.5) = +0.5
        """

        mapping = pd.DataFrame(
            {
                "batch": [
                    "B1",
                    "B1",
                    "B1",
                    "B1",
                    "B2",
                    "B2",
                    "B2",
                    "B2",
                ],
                "QC Lot": [
                    "L1",
                    "L1",
                    "L2",
                    "L2",
                    "L1",
                    "L1",
                    "L2",
                    "L2",
                ],
                "is_reference": [
                    True,
                    False,
                    True,
                    False,
                    True,
                    False,
                    True,
                    False,
                ],
            },
            index=[
                "S1",
                "S2",
                "S3",
                "S4",
                "S5",
                "S6",
                "S7",
                "S8",
            ],
        )

        row = pd.Series(
            [
                2.0,
                0.2,
                1.0,
                0.3,  # B1
                1.0,
                0.4,
                0.0,
                0.5,  # B2
            ],
            index=mapping.index,
        )

        expected = pd.Series(
            [
                1.5,
                -0.3,
                0.5,
                -0.2,  # B1 corrected (-0.5 shift)
                1.5,
                0.9,
                0.5,
                1.0,  # B2 corrected (+0.5 shift)
            ],
            index=row.index,
            name="intensity corrected",
        )

        result = bn.row_wise_normalize(row, mapping)

        pd.testing.assert_series_equal(result.round(6), expected.round(6))
