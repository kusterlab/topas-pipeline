from pathlib import Path

import pytest
import pandas as pd
import numpy as np

import topas_pipeline.clinical_tools as ct


class TestPhosphoAnnot:
    def test_annotates_phospho_sites_correctly(self, mocker):
        mock_df = pd.DataFrame(
            {
                "Modified sequence": ["seq1"],
                "PSP_LT_LIT": ["1;2"],
                "PSP_MS_LIT": ["3;4"],
                "PSP_MS_CST": ["5;6"],
                "PSP_MS_CST": ["5;6"],
            }
        )
        mock_add_positions = mocker.patch(
            "topas_pipeline.clinical_tools.pa.addPeptideAndPsitePositions", return_value=mock_df
        )
        mock_add_substrates = mocker.patch(
            "topas_pipeline.clinical_tools.pa.addPSPKinaseSubstrateAnnotations",
            return_value=mock_df,
        )
        mock_add_annotation = mocker.patch(
            "topas_pipeline.clinical_tools.pa.addPSPAnnotations", return_value=mock_df
        )
        mock_add_regulatory_annotations = mocker.patch(
            "topas_pipeline.clinical_tools.pa.addPSPRegulatoryAnnotations", return_value=mock_df
        )

        df = pd.DataFrame({"Modified sequence": ["seq1"]})
        extra_kinase_annot = False
        pspFastaFile = Path("pspFastaFile")
        pspKinaseSubstrateFile = Path("pspKinaseSubstrateFile")
        pspAnnotationFile = Path("pspAnnotationFile")
        pspRegulatoryFile = Path("pspRegulatoryFile")

        result_df = ct.phospho_annot(
            df,
            extra_kinase_annot,
            pspFastaFile,
            pspKinaseSubstrateFile,
            pspAnnotationFile,
            pspRegulatoryFile,
        )

        assert not result_df.empty
        assert mock_add_positions.call_count == 2
        mock_add_substrates.assert_called_once()
        mock_add_annotation.assert_called_once()
        mock_add_regulatory_annotations.assert_called_once()


class TestAddPspUrls:
    def test_add_psp_urls_correct_columns(self):
        data = {
            "Start positions": ["1", "2;4;2"],
            "Proteins": ["P12345", "P67890-2;P67890;P98765"],
        }
        df = pd.DataFrame(data)
        result = ct.add_psp_urls(df)

        expected_result = pd.DataFrame(
            {
                "Start positions": ["1", "2;4;2"],
                "Proteins": ["P12345", "P67890-2;P67890;P98765"],
                "PSP_URL": [
                    '=HYPERLINK("https://www.phosphosite.org/uniprotAccAction?id=P12345")',
                    '=HYPERLINK("https://www.phosphosite.org/uniprotAccAction?id=P67890")',
                ],
                "PSP_URL_extra": [
                    [],
                    [
                        "https://www.phosphosite.org/uniprotAccAction?id=P67890-2",
                        "https://www.phosphosite.org/uniprotAccAction?id=P98765",
                    ],
                ],
            }
        )
        pd.testing.assert_frame_equal(result, expected_result)


class TestAddUrlColumn:
    def test_valid_start_positions_and_proteins(self):
        row = ("1;2;3", "P12345;P67890;P54321")
        main_url, urls = ct.add_url_columns(row)
        assert (
            main_url
            == '=HYPERLINK("https://www.phosphosite.org/uniprotAccAction?id=P12345")'
        )
        assert urls == [
            "https://www.phosphosite.org/uniprotAccAction?id=P67890",
            "https://www.phosphosite.org/uniprotAccAction?id=P54321",
        ]


class TestProtBasketAnnotation:
    def test_annotates_dataframe_correctly(self, mocker):
        import pandas as pd
        from pathlib import Path
        from topas_pipeline.clinical_tools import prot_basket_annotation

        # Mocking the read_basket_annotation_generation functions
        mocker.patch(
            "topas_pipeline.clinical_tools.read_basket_annotation_generation1",
            return_value={"GeneA": {"basket": "Basket1", "weight": 0.5}},
        )
        mocker.patch(
            "topas_pipeline.clinical_tools.read_basket_annotation_generation2",
            return_value={"GeneA": {"basket": "Basket2", "weight": 0.7}},
        )
        mocker.patch(
            "topas_pipeline.clinical_tools.read_basket_annotation_generation4",
            return_value={"GeneA": {"basket": "Basket4", "weight": 0.9}},
        )

        # Creating a sample dataframe
        df = pd.DataFrame({"Gene names": ["GeneA", "GeneB"]})
        df = df.set_index("Gene names")

        # Calling the function under test
        result_df, annot_dict = prot_basket_annotation(
            df, "path/to/annotation/file", "fp", "basket"
        )

        expected_result_df = pd.DataFrame(
            {"basket": ["Basket1", ""], "basket_weights": [0.5, ""]},
            index=pd.Series(["GeneA", "GeneB"], name="Gene names"),
        )
        # Asserting the results
        assert "basket" in result_df.columns
        assert "basket_weights" in result_df.columns
        pd.testing.assert_frame_equal(result_df, expected_result_df)


class TestCreateIdentifierToBasketDict:
    def test_process_fp_data_type(self):
        data = {
            "gene": ["gene1", "gene2", "gene1"],
            "LEVEL": ["expression", "kinase activity", "expression"],
            "basket": ["basket1", "basket2", "basket3"],
        }
        df = pd.DataFrame(data)
        expected_output = {
            "gene1": {"LEVEL": "expression;expression", "basket": "basket1;basket3"},
            "gene2": {"LEVEL": "kinase activity", "basket": "basket2"},
        }

        result = ct.create_identifier_to_basket_dict(
            df, data_type="fp", identifier_type="gene"
        )
        assert result == expected_output

    def test_process_pp_data_type(self):
        data = {
            "gene": ["gene1", "gene2", "gene1"],
            "LEVEL": [
                "phosphorylation",
                "kinase activity",
                "important phosphorylation",
            ],
            "basket": ["basket1", "basket2", "basket3"],
        }
        df = pd.DataFrame(data)
        expected_output = {
            "gene1": {
                "LEVEL": "phosphorylation;important phosphorylation",
                "basket": "basket1;basket3",
            },
            "gene2": {"LEVEL": "kinase activity", "basket": "basket2"},
        }

        result = ct.create_identifier_to_basket_dict(
            df, data_type="pp", identifier_type="gene"
        )
        assert result == expected_output


class TestReadBasketAnnotationGeneration4:
    def test_reads_and_processes_valid_excel_file(self, mocker):
        mock_excel_data = pd.DataFrame(
            {
                "BASKET": ["basket1  ", "basket2"],
                "SUBBASKET": ["sub1", "sub2"],
                "LEVEL": ["expression  ", "kinase activity"],
                "WEIGHT": [0.5, 0.8],
                "GENE NAME": ["gene1", "gene2  "],
            }
        )
        mocker.patch("pandas.read_excel", return_value=mock_excel_data)

        result = ct.read_basket_annotation_generation4("dummy_path.xlsx", "fp")

        assert result == {
            "gene1": {
                "LEVEL": "expression",
                "basket": "basket1",
                "sub_basket": "sub1",
                "subbasket_level": "sub1 - expression",
                "weight": "0.5",
            },
            "gene2": {
                "LEVEL": "kinase activity",
                "basket": "basket2",
                "sub_basket": "sub2",
                "subbasket_level": "sub2 - kinase activity",
                "weight": "0.8",
            },
        }
