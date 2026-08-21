import sys
from pathlib import Path

import numpy as np
import pandas as pd


from topas_pipeline import sample_metadata


def get_metadata_df(metadata_file: Path):
    # not sure why, but these samples were not in the metadata
    missing_samples_df = pd.DataFrame(
        {
            "pat_H021-PHLFH6-M4": {"Batch_No": 36, "TMT_Channel": 3},
            "pat_H021-QC7SDC-T2": {"Batch_No": 38, "TMT_Channel": 4},
            "pat_H021-UQBN7H-T2": {"Batch_No": 36, "TMT_Channel": 4},
            "pat_I076-042-9T1": {"Batch_No": 36, "TMT_Channel": 5},
        }
    ).T
    missing_samples_df

    metadata_df = sample_metadata.load(metadata_file)
    metadata_df = metadata_df.set_index("Sample name")
    metadata_df.index = "pat_" + metadata_df.index
    metadata_df = pd.concat([metadata_df, missing_samples_df])

    metadata_df["Batch_Type"] = (
        metadata_df["Batch_No"]
        .astype(str)
        .map(lambda x: "Workflow_Test" if not x.isdigit() else "Sarcoma")
    )
    metadata_df["channel_name"] = (
        "Reporter intensity corrected "
        + metadata_df["TMT_Channel"].astype(str)
        + " "
        + metadata_df["Batch_Type"]
        + "_Batch"
        + metadata_df["Batch_No"].astype(str)
    )

    return metadata_df


def get_patient_channels_to_sample_id_dict(metadata_df: Path):
    patient_channels_to_sample_id_dict = {
        v: k for k, v in metadata_df["channel_name"].to_dict().items()
    }
    return patient_channels_to_sample_id_dict


def get_protein_log10_fold_change(
    results_folder: Path, patient_channels_to_sample_id_dict: dict[str, str]
):
    fp_df = pd.read_csv(
        results_folder / "preprocessed_fp2.csv",
        index_col="Gene names",
        usecols=lambda x: x == "Gene names" or x.startswith("Reporter intensity"),
    )

    fp_protein_median_df = fp_df.copy()
    fp_protein_median_df.columns = fp_protein_median_df.columns.map(
        patient_channels_to_sample_id_dict
    )
    fp_protein_median_df = fp_protein_median_df.loc[
        :, fp_protein_median_df.columns.notna()
    ]
    fp_protein_median_df = fp_protein_median_df.median(axis=1)

    fp_centered_df = fp_df.sub(fp_protein_median_df, axis=0)

    fp_centered_long_df = convert_wide_to_long(fp_centered_df)
    return fp_centered_long_df


def convert_wide_to_long(df: pd.DataFrame):
    idx = df.columns.str.rsplit(" ", n=1, expand=True)
    idx.names = ["TMT Channel", "Batch"]
    long_df = df.set_axis(idx, axis=1).stack()
    long_df = long_df.loc[
        :, [f"Reporter intensity corrected {i}" for i in range(1, 12)]
    ]
    return long_df


def load_evidence_df(results_folder: Path):
    # 4 minutes
    evidence_df = pd.read_csv(results_folder / "evidence.txt", sep="\t")
    return evidence_df


def get_rare_miscleaved_evidence_df(evidence_df: pd.DataFrame):
    evidence_df["Miscleavages"] = evidence_df["Modified sequence"].str[:-2].str.count("K|R")
    
    miscleaved_evidence_df = evidence_df[evidence_df["Miscleavages"] >= 1]
    miscleaved_evidence_df = miscleaved_evidence_df.groupby(
        ["Modified sequence", "Gene names", "Batch"]
    ).agg("sum")

    miscleaved_occurrence_evidence_df = (
        miscleaved_evidence_df.filter(like="Reporter intensity corrected")
        .replace(0, np.nan)
        .groupby("Modified sequence")
        .agg("count")
        .sum(axis=1)
    )

    miscleaved_evidence_df = miscleaved_evidence_df.join(
        miscleaved_occurrence_evidence_df.to_frame(name="occurrence")
    )

    tmt_columns = [f"Reporter intensity corrected {i}" for i in range(1, 12)]
    rare_miscleaved_evidence_df = miscleaved_evidence_df.loc[
        miscleaved_evidence_df["occurrence"] < 250, tmt_columns
    ]
    rare_miscleaved_evidence_df = np.log10(
        rare_miscleaved_evidence_df.replace(0, np.nan)
    )
    return rare_miscleaved_evidence_df


def normalize_by_fp(
    rare_miscleaved_evidence_df: pd.DataFrame, fp_centered_long_df: pd.DataFrame
):
    rare_miscleaved_evidence_fp_normalized_df = rare_miscleaved_evidence_df.sub(
        fp_centered_long_df, axis=0
    ).dropna(how="all")
    return rare_miscleaved_evidence_fp_normalized_df


def compute_miscleavage_index(rare_miscleaved_evidence_fp_normalized_df: pd.DataFrame):
    # Compute median intensity per sample
    rare_miscleaved_summed_intensity_df = (
        rare_miscleaved_evidence_fp_normalized_df.groupby("Batch").agg("median")
    )

    # normalize to ref channels
    rare_miscleaved_summed_intensity_df = rare_miscleaved_summed_intensity_df.subtract(
        rare_miscleaved_summed_intensity_df.loc[
            :, ["Reporter intensity corrected 10", "Reporter intensity corrected 11"]
        ].mean(axis=1),
        axis=0,
    )

    # convert to wide format
    rare_miscleaved_summed_intensity_wide_df = (
        rare_miscleaved_summed_intensity_df.unstack().to_frame(name="miscleavage_index")
    )
    rare_miscleaved_summed_intensity_wide_df.index = (
        rare_miscleaved_summed_intensity_wide_df.index.to_flat_index().str.join(" ")
    )

    return rare_miscleaved_summed_intensity_wide_df


def main(argv):
    results_folder = Path(argv[0])
    metadata_file = Path(argv[1])

    metadata_df = get_metadata_df(metadata_file)
    patient_channels_to_sample_id_dict = get_patient_channels_to_sample_id_dict(
        metadata_df
    )

    fp_centered_long_df = get_protein_log10_fold_change(
        results_folder, patient_channels_to_sample_id_dict
    )

    evidence_df = load_evidence_df(results_folder)
    rare_miscleaved_evidence_df = get_rare_miscleaved_evidence_df(evidence_df)

    rare_miscleaved_evidence_fp_normalized_df = normalize_by_fp(
        rare_miscleaved_evidence_df, fp_centered_long_df
    )

    rare_miscleaved_summed_intensity_wide_df = compute_miscleavage_index(
        rare_miscleaved_evidence_fp_normalized_df
    )
    rare_miscleaved_summed_intensity_wide_df = (
        rare_miscleaved_summed_intensity_wide_df.merge(
            metadata_df[["channel_name", "Program"]],
            left_index=True,
            right_on="channel_name",
        )
    )
    rare_miscleaved_summed_intensity_wide_df.to_csv(
        results_folder / "miscleavage_index.tsv", sep="\t"
    )


if __name__ == "__main__":
    main(sys.argv[1:])