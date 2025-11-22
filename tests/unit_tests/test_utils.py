import pandas as pd
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


def test_merge_by_delimited_field_basic():
    # Input dataframe
    df = pd.DataFrame(
        {"genes": ["TP53;BRCA1", "EGFR;KRAS;BRAF"], "other_annot": ["AAA", "BBB"]}
    )

    # Lookup dataframe
    mapping = pd.DataFrame(
        {
            "genes": ["TP53", "BRCA1", "EGFR", "KRAS", "BRAF"],
            "pathway": ["P53", "BRCA", "EGFR", "MAPK", "RAF"],
        }
    )

    # Expected outcome: each gene maps to its pathway, merged back by id
    expected = pd.DataFrame(
        {
            "genes": ["TP53;BRCA1", "EGFR;KRAS;BRAF"],
            "pathway": ["P53;BRCA", "EGFR;MAPK;RAF"],
            "other_annot": ["AAA", "BBB"],
        }
    )

    result = utils.merge_by_delimited_field(df, mapping, field_name="genes")

    # Check that results match
    pd.testing.assert_frame_equal(
        result, expected, check_like=True, check_index_type=False
    )


def test_merge_by_delimited_field_index_col():
    # Input dataframe
    df = pd.DataFrame(
        {"genes": ["TP53;BRCA1", "EGFR;KRAS;BRAF"], "other_annot": ["AAA", "BBB"]}
    ).set_index("genes")

    # Lookup dataframe
    mapping = pd.DataFrame(
        {
            "genes": ["TP53", "BRCA1", "EGFR", "KRAS", "BRAF"],
            "pathway": ["P53", "BRCA", "EGFR", "MAPK", "RAF"],
        }
    )

    # Expected outcome: each gene maps to its pathway, merged back by id
    expected = pd.DataFrame(
        {
            "genes": ["TP53;BRCA1", "EGFR;KRAS;BRAF"],
            "pathway": ["P53;BRCA", "EGFR;MAPK;RAF"],
            "other_annot": ["AAA", "BBB"],
        }
    )

    result = utils.merge_by_delimited_field(df, mapping, field_name="genes")

    # Check that results match
    pd.testing.assert_frame_equal(
        result, expected, check_like=True, check_index_type=False
    )


def test_merge_by_delimited_field_inplace():
    # Input dataframe
    df = pd.DataFrame(
        {"genes": ["TP53;BRCA1", "EGFR;KRAS;BRAF"], "other_annot": ["AAA", "BBB"]}
    )

    # Lookup dataframe
    mapping = pd.DataFrame(
        {
            "genes": ["TP53", "BRCA1", "EGFR", "KRAS", "BRAF"],
            "pathway": ["P53", "BRCA", "EGFR", "MAPK", "RAF"],
        }
    )

    # Expected outcome: each gene maps to its pathway, merged back by id
    expected = pd.DataFrame(
        {
            "genes": ["TP53;BRCA1", "EGFR;KRAS;BRAF"],
            "pathway": ["P53;BRCA", "EGFR;MAPK;RAF"],
            "other_annot": ["AAA", "BBB"],
        }
    )

    utils.merge_by_delimited_field(df, mapping, field_name="genes", inplace=True)

    # Check that results match
    pd.testing.assert_frame_equal(df, expected, check_like=True, check_index_type=False)


def test_merge_by_delimited_field_duplicates():
    # Input dataframe
    df = pd.DataFrame(
        {
            "genes": ["TP53;BRCA1", "TP53;BRCA1", "EGFR;KRAS;BRAF"],
            "other_annot": ["AAA", "BBB", "CCC"],
        }
    )

    # Lookup dataframe
    mapping = pd.DataFrame(
        {
            "genes": ["TP53", "BRCA1", "EGFR", "KRAS", "BRAF"],
            "pathway": ["P53", "BRCA", "EGFR", "MAPK", "RAF"],
        }
    )

    # Expected outcome: each gene maps to its pathway, merged back by id
    expected = pd.DataFrame(
        {
            "genes": ["TP53;BRCA1", "TP53;BRCA1", "EGFR;KRAS;BRAF"],
            "pathway": ["P53;BRCA", "P53;BRCA", "EGFR;MAPK;RAF"],
            "other_annot": ["AAA", "BBB", "CCC"],
        }
    )

    result = utils.merge_by_delimited_field(df, mapping, field_name="genes")

    # Check that results match
    pd.testing.assert_frame_equal(
        result, expected, check_like=True, check_index_type=False
    )


def test_merge_by_delimited_field_handles_missing_values():
    df = pd.DataFrame({"genes": ["TP53;UNKNOWN"]})

    mapping = pd.DataFrame({"genes": ["TP53"], "pathway": ["P53"]})

    result = utils.merge_by_delimited_field(df, mapping, field_name="genes")

    # Missing entries should result in NaN in merged columns
    assert "UNKNOWN" in result["genes"].iloc[0]
    assert result["pathway"].iloc[0] == "P53"


def test_merge_by_delimited_field_empty_df():
    df = pd.DataFrame(columns=["id", "genes"]).set_index("id")
    mapping = pd.DataFrame(columns=["genes", "pathway"])

    result = utils.merge_by_delimited_field(df, mapping, "genes")

    assert isinstance(result, pd.DataFrame)
    assert result.empty