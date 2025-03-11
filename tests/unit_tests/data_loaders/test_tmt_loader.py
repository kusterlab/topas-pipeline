import pandas as pd
import numpy as np

import bin.data_loaders.tmt_loader as tmt_loader


def test_scale_ms1_to_reference():
    # Test data setup
    data = {
        "Modified sequence": ["pep1", "pep1", "pep2", "pep2"],
        "Charge": [2, 2, 3, 3],
        "Batch": ["01", "02", "01", "02"],
        "MS1": [100, 0, 300, 200],  # MS1 with a zero to test replacement with NaN
        "Reporter intensity corrected 1": [200, 700, 400, 300],
        "Reporter intensity corrected 2": [300, 300, 400, 500],
    }

    # pep1:
    # ref channel contribution: Batch1: 300/500=0.6; Batch 2: 300/1000=0.3
    # ms1 contributions ref channels: Batch1: 100*0.6=60; Batch 2: NaN
    # mean ms1 contributions ref channels: 60
    # ms1 after scaling to reference: Batch1: 60/0.6=100; Batch 2: 60/0.3=200

    # pep2:
    # ref channel contribution: Batch1: 400/800=0.5; Batch 2: 500/800=0.625
    # ms1 contributions ref channels: Batch1: 300*0.5=150; Batch 2: 200*0.625=125
    # mean ms1 contributions ref channels: 137.5
    # ms1 after scaling to reference: Batch1: 137.5/0.5=275; Batch 2: 137.5/0.625=220

    sample_annot = {
        "Sample name": ["AAA", "BBB", "CCC", "DDD"],
        "Batch Name": ["01", "01", "02", "02"],
        "TMT Channel": [1, 2, 1, 2],
        "is_reference": [False, True, False, True],
    }
    df = pd.DataFrame(data)
    sample_annotation_df = pd.DataFrame(sample_annot)
    ref_channels_df = sample_annotation_df.loc[sample_annotation_df["is_reference"], :]

    result_df = tmt_loader._scale_ms1_to_reference(df, ref_channels_df)

    # Assertions
    assert "MS1 corrected reference intensity" in result_df.columns
    assert "Mean MS1 corrected reference intensity" in result_df.columns

    expected_ms1_corrected = pd.Series([100.0, 200.0, 275.0, 220.0], name="MS1")
    pd.testing.assert_series_equal(result_df["MS1"], expected_ms1_corrected)


def test_scale_ms1_to_reference_and_redistribute():
    # Test data setup
    data = {
        "Modified sequence": ["pep1", "pep1", "pep2", "pep2"],
        "Charge": [2, 2, 3, 3],
        "Batch": ["01", "02", "01", "02"],
        "MS1": [100, 0, 300, 200],  # MS1 with a zero to test replacement with NaN
        "Reporter intensity corrected 1": [200, 700, 400, 300],
        "Reporter intensity corrected 2": [300, 300, 400, 500],
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

    result_df = tmt_loader._scale_ms1_to_reference(df, ref_channels_df)
    result_df = tmt_loader._redistribute_ms1_intensity(result_df)

    # check that the reference channels intensities for each peptide are equal after redistribution
    expected_df = pd.DataFrame(
        data={
            "Modified sequence": ["pep1", "pep1", "pep2", "pep2"],
            "Charge": [2, 2, 3, 3],
            "Batch": ["01", "02", "01", "02"],
            "MS1": [100.0, 200.0, 275.0, 220.0],
            "Reporter intensity corrected 1": [40.0, 140.0, 137.5, 82.5],
            "Reporter intensity corrected 2": [60.0, 60.0, 137.5, 137.5],
            "MS1 corrected reference intensity": [60.0, np.nan, 150.0, 125.0],
            "Mean MS1 corrected reference intensity": [60, 60, 137.5, 137.5],
        }
    )
    pd.testing.assert_frame_equal(result_df, expected_df)


def test_impute_ms1_intensity():
    # Test data setup
    data = {
        "Modified sequence": ["pep1", "pep1", "pep2", "pep2"],
        "Charge": [2, 2, 3, 3],
        "Batch": ["01", "02", "01", "02"],
        "MS1": [100, 0, 300, 200],  # MS1 with a zero to test replacement with NaN
        "Reporter intensity corrected 1": [200, 700, 400, 300],
        "Reporter intensity corrected 2": [300, 300, 400, 500],
    }

    # pep1:
    # ref channel contribution: Batch1: 300/500=0.6; Batch 2: 300/1000=0.3
    # ms1 contributions ref channels: Batch1: 100*0.6=60; Batch 2: NaN
    # mean ms1 contributions ref channels: 60
    # ms1 after scaling to reference: Batch1: 60/0.6=100; Batch 2: 60/0.3=200

    # pep2:
    # ref channel contribution: Batch1: 400/800=0.5; Batch 2: 500/800=0.625
    # ms1 contributions ref channels: Batch1: 300*0.5=150; Batch 2: 200*0.625=125
    # mean ms1 contributions ref channels: 137.5
    # ms1 after scaling to reference: Batch1: 137.5/0.5=275; Batch 2: 137.5/0.625=220

    sample_annot = {
        "Sample name": ["AAA", "BBB", "CCC", "DDD"],
        "Batch Name": ["01", "01", "02", "02"],
        "TMT Channel": [1, 2, 1, 2],
        "is_reference": [False, True, False, True],
    }
    df = pd.DataFrame(data)
    sample_annotation_df = pd.DataFrame(sample_annot)
    ref_channels_df = sample_annotation_df.loc[sample_annotation_df["is_reference"], :]

    result_df = tmt_loader._impute_ms1_intensity(df, ref_channels_df)

    # Assertions
    assert "MS1 corrected reference intensity" in result_df.columns
    assert "Mean MS1 corrected reference intensity" in result_df.columns

    expected_ms1_corrected = pd.Series([100.0, 200.0, 300.0, 200.0], name="MS1")
    pd.testing.assert_series_equal(result_df["MS1"], expected_ms1_corrected)