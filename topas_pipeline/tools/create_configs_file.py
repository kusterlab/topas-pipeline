"""

Description:
Script to create input sample annotation file for pipeline from metadata file.
The location of MQ searches are given by Cohort column info and is either hard coded to be Sarcoma (for patients) or Workflow_Test (for models) based on Batch name.
Reference channels are added for each batch as TMT channel 10 and 11.

Usage:

    python tools/create_sample_annotation.py Metadata_file.xlsx output_folder

"""

import sys
import re
import argparse
from pathlib import Path

import pandas as pd

from topas_pipeline import config


def main(argv):
    parser = argparse.ArgumentParser(
        description="Script to create input sample annotation file for pipeline from metadata file."
    )
    parser.add_argument(
        "--config-file-template",
        required=True,
        help="Config file to use as template. The result_folder, sample_annotation and metadata_annotation fields will be overwritten.",
    )
    parser.add_argument(
        "--metadata-file", required=True, help="Path to the metadata Excel file"
    )
    parser.add_argument(
        "--results-folder",
        required=True,
        help="Path to the folder where to write the pipeline results, should be of the format 2025.11.04_xxx",
    )
    parser.add_argument(
        "--qc-lot-mapping-file",
        help="Path to Excel file with two mandatory columns: 'Batch Name' and 'QC Lot'",
    )

    args = parser.parse_args(argv)

    metadata_file = args.metadata_file
    qc_lot_mapping_file = args.qc_lot_mapping_file
    results_folder = Path(args.results_folder)
    output_file_suffix = results_folder.name

    results_folder.mkdir(exist_ok=True)

    sample_annotation_file = create_sample_annotation(
        metadata_file,
        results_folder,
        qc_lot_mapping_file,
        output_file_suffix,
    )

    configs = config.load(args.config_file_template)
    configs.results_folder = str(results_folder)
    configs.sample_annotation = sample_annotation_file
    configs.metadata_annotation = metadata_file

    with open(results_folder / "configs.toml", "w") as toml_file:
        toml_file.write(configs.astoml())


def create_sample_annotation(
    metadata_file: str,
    output_folder: str,
    qc_lot_mapping_file: str,
    output_file_suffix: str,
):
    # Rest of your logic goes here
    print(f"Metadata file: {metadata_file}")
    print(f"Output folder: {output_folder}")
    print(f"QC lot mapping file: {qc_lot_mapping_file}")

    metadata_df = pd.read_excel(metadata_file, engine="openpyxl")
    metadata_df = metadata_df[metadata_df["QC"] != "exclude"]
    metadata_df = metadata_df.rename(
        columns={"Batch_No": "Batch Name", "TMT_Channel": "TMT Channel"}
    )

    new_ref_rows = []
    # Iterate through each unique Batch Name
    for batch in metadata_df["Batch Name"].unique():
        for channel in [10, 11]:
            # Create a new row as a dictionary
            new_row = {
                "Batch Name": batch,
                "Sample name": f"ref_channel_{channel}_batch{batch}",
                "TMT Channel": channel,
                "QC": "passed",
                "is_reference": True,
            }
            # Append the new row to the list
            new_ref_rows.append(new_row)

    # Convert the list of new rows into a DataFrame
    ref_channel_df = pd.DataFrame(new_ref_rows)

    if qc_lot_mapping_file:
        qc_lot_mapping_df = pd.read_excel(qc_lot_mapping_file, dtype=str)

        qc_lot_mapping_df = pd.DataFrame(
            [x for _, row in qc_lot_mapping_df.iterrows() for x in expand_rows(row)]
        )
        duplicates = qc_lot_mapping_df[
            qc_lot_mapping_df.duplicated(
                subset=["Batch Name", "TMT Channel"], keep=False
            )
        ]
        if len(duplicates) > 0:
            raise ValueError(
                f"Found duplicate entries in QC lot mapping file:\n{duplicates}"
            )

        ref_channel_df = ref_channel_df.merge(
            qc_lot_mapping_df,
            on=["Batch Name", "TMT Channel"],
            how="left",
            validate="1:1",
        )

        if ref_channel_df["QC Lot"].hasnans:
            raise ValueError(
                f"Some QC channels did not have QC Lots:\n{ref_channel_df[ref_channel_df['QC Lot'].isna()]}"
            )

    # Concatenate the new rows DataFrame with the original DataFrame
    metadata_df = pd.concat([metadata_df, ref_channel_df], ignore_index=True)

    metadata_df["Cohort"] = "Sarcoma"

    def update_cohort(row):
        batch_name = str(row["Batch Name"])  # Convert the value to string
        # Check if the batch name contains both letters and numbers
        if re.search(r"[A-Za-z]", batch_name) and re.search(r"\d", batch_name):
            return "Workflow_Test"
        return row["Cohort"]

    metadata_df["Cohort"] = metadata_df.apply(update_cohort, axis=1)

    metadata_df = metadata_df.sort_values(["Batch Name", "TMT Channel"])

    output_file = f"{output_folder}/sample_annotation_{output_file_suffix}.csv"
    metadata_df.to_csv(output_file, index=False)

    print(f"Written sample annotation file to: {output_file}")

    return output_file


def expand_rows(row: pd.Series):
    batch = row["Batch Name"]
    qc = row["QC Lot"]

    # Expand batch range or leave as string
    if "-" in batch:
        start, end = batch.split("-")
        batches = list(range(int(start), int(end) + 1))
    else:
        batches = [batch]  # keep string like 'CL15'

    # If mixed, create per-channel rows
    mixed_match = re.findall(r"channel (\d+)=(\d+)", qc)
    for b in batches:
        if mixed_match:
            for channel, lot in mixed_match:
                yield {"Batch Name": b, "TMT Channel": int(channel), "QC Lot": int(lot)}
        else:
            # Same QC for all channels 10 and 11
            for ch in [10, 11]:
                yield {"Batch Name": b, "TMT Channel": ch, "QC Lot": int(qc)}


if __name__ == "__main__":
    main(sys.argv[1:])
