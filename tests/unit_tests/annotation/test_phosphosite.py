from pathlib import Path

import pandas as pd
import numpy as np

from topas_pipeline.annotation import phosphosite
from topas_pipeline import config


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
            "psite_annotation.addPeptideAndPsitePositions",
            return_value=mock_df,
        )
        mock_add_substrates = mocker.patch(
            "psite_annotation.addPSPKinaseSubstrateAnnotations",
            return_value=mock_df,
        )
        mock_add_annotation = mocker.patch(
            "psite_annotation.addPSPAnnotations", return_value=mock_df
        )
        mock_add_regulatory_annotations = mocker.patch(
            "psite_annotation.addPSPRegulatoryAnnotations",
            return_value=mock_df,
        )

        df = pd.DataFrame({"Modified sequence": ["seq1"]})
        clinic_proc_config = config.ClinicProc(
            extra_kinase_annot="",
            pspFastaFile=Path("pspFastaFile"),
            pspKinaseSubstrateFile=Path("pspKinaseSubstrateFile"),
            pspAnnotationFile=Path("pspAnnotationFile"),
            pspRegulatoryFile=Path("pspRegulatoryFile"),
            prot_baskets="",
        )

        result_df = phosphosite.add_phospho_annotations(
            df, clinic_proc_config=clinic_proc_config
        )

        assert not result_df.empty
        assert mock_add_positions.call_count == 2
        assert mock_add_substrates.call_count == 2
        mock_add_annotation.assert_called_once()
        mock_add_regulatory_annotations.assert_called_once()
