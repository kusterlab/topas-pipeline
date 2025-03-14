import pandas as pd
import numpy as np
import pytest

import topas_pipeline.utils as utils


def test_get_ref_channels():
    # Test data setup
    data = {
        "Modified sequence": ["pep1", "pep1", "pep2", "pep2"],
        "Batch": ["01", "02", "01", "02"],
        "MS1": [100.0, 0.0, 300.0, 200.0],
        "Reporter intensity corrected 1": [200.0, 700.0, 400.0, 300.0],
        "Reporter intensity corrected 2": [300.0, 300.0, 400.0, 500.0],
    }
    sample_annot = {
        "Sample name": ["AAA", "BBB", "CCC", "DDD"],
        "Batch Name": ["01", "01", "02", "02"],
        "TMT Channel": [1, 2, 1, 2],
        "is_reference": [False, True, False, True],
    }
    df = pd.DataFrame(data)
    sample_annotation_df = pd.DataFrame(sample_annot)
    ref_channels_df = sample_annotation_df.loc[sample_annotation_df["is_reference"], :]

    result_df = utils.get_ref_channels(df, ref_channels_df)

    # check that only the reference channel intensities are retained and the rest are NaNs
    expected_df = pd.DataFrame(
        data={
            "Reporter intensity corrected 2": [300.0, 300.0, 400.0, 500.0],
        }
    )
    pd.testing.assert_frame_equal(result_df, expected_df)


def test_get_ref_channels():
    # Test data setup
    data = {
        "Modified sequence": ["pep1", "pep1", "pep2", "pep2"],
        "Batch": ["01", "02", "01", "02"],
        "MS1": [100.0, 0.0, 300.0, 200.0],
        "Reporter intensity corrected 1": [200.0, 700.0, 400.0, 300.0],
        "Reporter intensity corrected 2": [300.0, 300.0, 400.0, 500.0],
    }
    sample_annot = {
        "Sample name": ["AAA", "BBB", "CCC", "DDD"],
        "Batch Name": ["01", "01", "02", "02"],
        "TMT Channel": [1, 2, 1, 2],
        "is_reference": [True, False, False, True],
    }
    df = pd.DataFrame(data)
    sample_annotation_df = pd.DataFrame(sample_annot)
    ref_channels_df = sample_annotation_df.loc[sample_annotation_df["is_reference"], :]

    with pytest.raises(
        ValueError,
        match="Datasets where reference channels differ between batches are not yet supported.",
    ):
        utils.get_ref_channels(df, ref_channels_df)
