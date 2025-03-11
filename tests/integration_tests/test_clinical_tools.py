import pandas as pd
import numpy as np
from pathlib import Path

from bin import clinical_tools, config

CONFIG_FILE_PATH = './data/test_config.json'


def test_phospho_annot():
    """
    This only works after replacing all empty PSP files in location /data/annotation_files
    Further after first succesful run of pipeline add path of results_folder to test_config.json
    And last add the annotated preprocessed_pp.csv file (annot_pp.csv) to integration_tests folder
    """
    configs = config.load(CONFIG_FILE_PATH)

    preprocessed_pp_file = Path(configs["results_folder"]) / "preprocessed_pp.csv"
    reference_result_file = Path(configs["results_folder"]) / "integration_tests" / "preprocessed_pp_phospho_annot.csv"

    index_col = "Modified sequence"
    keep_default_na = False
    
    preprocessed_pp = pd.read_csv(
        preprocessed_pp_file, index_col=index_col, keep_default_na=keep_default_na
    )

    patient_cols = preprocessed_pp.filter(regex=r"^pat_").columns
    preprocessed_pp[patient_cols] = (preprocessed_pp[patient_cols].replace("", np.nan).astype("float"))

    preprocessed_pp = clinical_tools.phospho_annot(
        preprocessed_pp,
        extra_kinase_annot=configs["extra_kinase_annot"],
        pspFastaFile=configs["clinic_proc"]["pspFastaFile"],
        pspKinaseSubstrateFile=configs["clinic_proc"]["pspKinaseSubstrateFile"],
        pspAnnotationFile=configs["clinic_proc"]["pspAnnotationFile"],
        pspRegulatoryFile=configs["clinic_proc"]["pspRegulatoryFile"],
    )

    # uncomment this to create new reference file
    preprocessed_pp.to_csv(reference_result_file, sep="\t", index=False)

    preprocessed_pp_reference = pd.read_csv(
        reference_result_file,
        sep="\t",
        dtype={"PSP_LT_LIT": "object", "PSP_MS_LIT": "object", "PSP_MS_CST": "object"},
    )
    non_patient_cols = [
        c for c in preprocessed_pp_reference.columns if c not in patient_cols
    ]
    preprocessed_pp_reference[non_patient_cols] = preprocessed_pp_reference[
        non_patient_cols
    ].fillna("")
    
    pd.testing.assert_frame_equal(
        preprocessed_pp,
        preprocessed_pp_reference,
        check_like=True,
        check_dtype=False,
        check_exact=False,
    )


if __name__ == "__main__":
    test_phospho_annot()
