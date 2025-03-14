import pandas as pd
from unittest.mock import patch
from topas_pipeline.picked_group import remap_gene_names


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
