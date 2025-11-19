from unittest.mock import MagicMock

import pandas as pd
import pytest
import numpy as np

from picked_group_fdr.results import ProteinGroupResult, ProteinGroupResults
from picked_group_fdr.protein_groups import ProteinGroups
from picked_group_fdr.precursor_quant import PrecursorQuant

from topas_pipeline.preprocess import picked_group

# from unittest.mock import patch
#
# @patch("picked_group_fdr.digest.getProteinAnnotations")
# def test_remap_gene_names(mock_getProteinAnnotations):
#     # Mock data
#     fasta_file = "test.fasta"
#     df = pd.DataFrame(
#         {
#             "Proteins": ["A;B;C", "D;E;F", "G", "H", ""],
#             "Gene names": [
#                 "GeneOldA;GeneOldB;GeneOldC",
#                 "GeneOldD;GeneOldE;GeneOldF",
#                 "GeneOldG",
#                 "GeneOldH",
#                 "GeneOldI",
#             ],
#         }
#     )

#     # Mock the return values for the mocked functions
#     mock_getProteinAnnotations.return_value = {
#         "A": ("A", "GeneA", "HeaderA"),
#         "B": ("B", "GeneB", "HeaderB"),
#         "C": ("C", "GeneC", "HeaderC"),
#         "D": ("D", "", "HeaderD"),
#         "E": ("E", "GeneE", "HeaderE"),
#         "F": ("F", "GeneF", "HeaderF"),
#         "G": ("G", "GeneG", "HeaderG"),
#     }

#     # Call the function
#     result = remap_gene_names(df, fasta_file)

#     # Check the output
#     expected_df = pd.DataFrame(
#         {
#             "Proteins": ["A;B;C", "D;E;F", "G", "H", ""],
#             "Gene names": [
#                 "GeneA;GeneB;GeneC",
#                 "GeneE;GeneF",
#                 "GeneG",
#                 "GeneOldH",
#                 "GeneOldI",
#             ],
#         }
#     )
#     pd.testing.assert_frame_equal(result, expected_df)


@pytest.fixture
def sample_dataframe():
    """Creates a sample dataframe similar to expected input."""
    data = {
        "Modified sequence": ["_PEPTIDE1_", "_PEPTIDE2_", "_PEPTIDE3_"],
        "Charge": [2, 3, 2],
        "Batch": ["A", "A", "B"],
        "Fraction": [1, 2, 1],
        "PEP": [0.01, 0.02, 0.03],
        "id": [101, 102, 103],
        "Reporter intensity corrected 1": [100.0, 200.0, np.nan],
        "Reporter intensity corrected 2": [50.0, 0.0, 300.0],  # One zero intensity
        "Reporter intensity corrected 3": [50.0, np.nan, 300.0],  # One zero intensity
    }
    return pd.DataFrame(data)


def test_add_peptides_to_protein_groups(sample_dataframe):
    """Tests if the function correctly adds peptides to protein groups."""
    peptideToProteinMap = {
        "PEPTIDE1": {"P12345"},
        "PEPTIDE2": {"P67890"},
        "PEPTIDE3": set(),
    }
    proteinGroups = ProteinGroups([["P12345"], ["P67890"]])
    proteinGroups.create_index()
    proteinGroupResults = ProteinGroupResults([ProteinGroupResult("P12345"), ProteinGroupResult("P67890")])

    # Call function
    result = picked_group.add_peptides_to_protein_groups(
        sample_dataframe, peptideToProteinMap, proteinGroupResults, proteinGroups
    )

    # Check if precursorQuants were added correctly
    assert (
        len(list(proteinGroupResults[0].precursorQuants)) == 3
    )  # Three valid intensities

    # Check if precursorQuants were added correctly
    assert (
        len(list(proteinGroupResults[1].precursorQuants)) == 1
    )  # One valid intensities (excluding zero and NaN)

    assert "Reporter intensity corrected 1 A" in result.experiments
    assert "Reporter intensity corrected 2 A" in result.experiments
