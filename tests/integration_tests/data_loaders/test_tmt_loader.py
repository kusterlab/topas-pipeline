from pathlib import Path

import numpy as np
import pandas as pd

from topas_pipeline.data_loaders.tmt_loader import TMTLoader
from topas_pipeline.preprocess.preprocess_tools import MQ_EVIDENCE_COLUMNS


def test_tmt_data_loader():
    evidence_files = [
        str(Path(__file__).parent / "Batch01_FP_CPTAC_UCEC" / "evidence.txt")
    ]

    # # generate minimum input file
    # evidence_df = pd.read_csv(evidence_files[0], sep="\t", usecols=MQ_EVIDENCE_COLUMNS, nrows=10)
    # evidence_df.to_csv(evidence_files[0], sep="\t", index=False)

    sample_annot = {
        "Sample name": ["AAA", "BBB", "CCC", "DDD"],
        "Batch": ["data_loaders_Batch01", "data_loaders_Batch01", "data_loaders_Batch01", "data_loaders_Batch01"],
        "TMT Channel": [1, 2, 3, 4],
        "is_reference": [False, False, False, True],
        "QC Lot": [np.nan, np.nan, np.nan, 1.0],
    }
    sample_annotation_df = pd.DataFrame(sample_annot)
    ref_channels_df = sample_annotation_df.loc[sample_annotation_df["is_reference"], :]

    data_loader = TMTLoader(evidence_files)
    all_batches = data_loader.load_data(MQ_EVIDENCE_COLUMNS)
    df = pd.concat(all_batches, ignore_index=True)

    df, correction_factors = data_loader.median_centering_within_batch(df)

    df, correction_factors = data_loader.median_centering_ms1(df)

    df = data_loader.impute_ms1_intensity(df, ref_channels_df)

    df = data_loader.redistribute_ms1_intensity(df)

    print(df)


if __name__ == "__main__":
    test_tmt_data_loader()
