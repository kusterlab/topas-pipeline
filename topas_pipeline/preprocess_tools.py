from builtins import ValueError
import re
import os
from glob import glob
from typing import List, Dict, Tuple, Union, Callable
from pathlib import Path
import logging

import pandas as pd
import numpy as np

from . import utils
from . import sample_metadata
from . import sample_annotation
from . import identification_metadata as id_meta

logger = logging.getLogger(__name__)
MQ_EVIDENCE_COLUMNS = [
    "Modifications",
    "Modified sequence",
    "Proteins",
    "Leading proteins",
    "Gene names",
    "Type",
    "Raw file",
    "Fraction",
    "Experiment",
    "Charge",
    "PEP",
    "MS/MS scan number",
    "Score",
    "Intensity",
    "Reverse",
    "Potential contaminant",
    "id",
] + [f"Reporter intensity corrected {i}" for i in range(1, 12)]

MQ_EVIDENCE_COLUMNS_TYPES = {
    "Modifications": "category",
    "Modified sequence": "",
    "Proteins": "",
    "Leading proteins": "",
    "Gene names": "",
    "Type": "",
    "Raw file": "",
    "Fraction": "",
    "Experiment": "category",
    "Charge": "",
    "PEP": "",
    "MS/MS scan number": "",
    "Score": "",
    "Intensity": "",
    "Reverse": "",
    "Potential contaminant": "",
    "id": "",
}


# TODO: move this to config validation
def check_annot(
    sample_annotation_file: str, metadata_annotation_file: str, in_metadata: Callable
):
    """Performs consistency checks of the sample annotation and metadata files.

    Args:
        sample_annotation_file: _description_
        metadata_annotation_file: _description_
        in_metadata: _description_

    Raises:
        ValueError: _description_
        ValueError: _description_
        ValueError: _description_
        ValueError: _description_

    Returns:
        _description_

    """
    sample_annot_df = sample_annotation.load_sample_annotation(sample_annotation_file)
    sample_annot_df_filtered = sample_annotation.filter_sample_annotation(
        sample_annot_df, remove_qc_failed=True, remove_replicates=True
    )
    sample_annot_df_filtered = sample_annot_df_filtered.reset_index()  # drop=True
    metadata_df = sample_metadata.load(metadata_annotation_file)

    # Check that there are no duplicates in neither patient annot and metadata
    if sample_annot_df_filtered["Sample name"].duplicated().any():
        duplicated = sample_annot_df_filtered[
            sample_annot_df_filtered["Sample name"].duplicated()
        ]
        logger.info(f"Duplicated sample(s) in sample annotation: {duplicated}")
        raise ValueError(f"Duplicated sample(s) in sample annotation: {duplicated}")

    elif metadata_df["Sample name"].duplicated().any():
        duplicates_samples = metadata_df["Sample name"][
            metadata_df["Sample name"].duplicated()
        ]
        duplicates_samples = "_".join(list(duplicates_samples))
        logger.info(
            f"The Duplicated sample(s): {duplicates_samples} in metadata annotation should be removed"
        )
        raise ValueError(
            f"The Duplicated sample(s): {duplicates_samples} in metadata annotation should be removed"
        )

    # TODO : needs to be refactored

    # Check for duplicates in batch, tmt_channel
    elif (
        sample_annot_df_filtered[["Cohort", "Batch Name", "TMT Channel"]]
        .duplicated()
        .any()
    ):
        duplicated = sample_annot_df_filtered[
            sample_annot_df_filtered[
                ["Cohort", "Batch Name", "TMT Channel"]
            ].duplicated()
        ]
        logger.info(
            f"Duplicated cohort, batch and tmt_channel in sample annotation: {duplicated}"
        )
        raise ValueError(
            f"Duplicated cohort, batch and tmt_channel in sample annotation: {duplicated}"
        )

    # Check if all given samples in patient annot is also in metadata
    # if sample_annot_df_filtered['Sample name'][~sample_annot_df_filtered['Sample name'].apply(lambda x: in_metadata(x)).isin(metadata_df['Sample name'])].any():
    #    missing_samples = sample_annot_df_filtered['Sample name'][~sample_annot_df_filtered['Sample name'].apply(lambda x: in_metadata(x)).isin(metadata_df['Sample name'])].values
    # raise ValueError(f'Sample(s) not found in metadata: {missing_samples}')
    #    print(f'Sample(s) not found in metadata: {missing_samples}')

    return sample_annot_df_filtered


def in_metadata(sample: str = "CHD"):
    # TODO: refactor this with a more generic function removing replicate suffixes
    if sample.split("-")[-1] == "R2":
        return "-".join(sample.split("-")[:-1])
    return sample


def remove_ref_empty_batch(
    df_with_ref: pd.DataFrame, sample_annotation_df: pd.DataFrame
):
    # remove ref for empty batches
    sample_names_drop = []
    for col in df_with_ref.columns:
        if "Reporter intensity corrected" in col:
            # TODO: handle case where "Batch" is not contained in column name
            # TODO: check if this works for batch numbers with 3 or more digits
            batch = re.search(r"Batch\d{1,2}", col).group()
            batch = re.findall(r"\d+", batch)[0]
            qc_passed = sample_annotation_df.loc[
                sample_annotation_df["Batch Name"] == int(batch), "QC"
            ]
            if "passed" not in qc_passed.values and "shaky" not in qc_passed.values:
                sample_names_drop.append(col)
    # samples have already been dropped but we want to drop the reporter intensity ones for this batch
    df_with_ref = df_with_ref.drop(columns=sample_names_drop)
    return df_with_ref


def get_files_by_type(
    raw_data_location: str,
    data_type: str,
    sample_annotation_df: pd.DataFrame,
    file_type: str,
):
    logger.info(f"Getting paths to {file_type} files for {data_type}")

    batches = sample_annotation.get_unique_batches(sample_annotation_df)
    found_files = get_data_location(raw_data_location, data_type, file_type=file_type)
    found_files = filter_evidence_files(found_files, data_type.upper(), batches)
    return found_files


def get_summary_files(raw_data_location, data_type, sample_annotation_file):
    sample_annotation_df = sample_annotation.load_sample_annotation(
        sample_annotation_file
    )
    return get_files_by_type(
        raw_data_location, data_type, sample_annotation_df, file_type="summary.txt"
    )


def get_evidence_files(
    sample_annotation_df: pd.DataFrame,
    raw_data_location: Union[str, Path],
    data_type: str,
) -> List:
    return get_files_by_type(
        raw_data_location, data_type, sample_annotation_df, file_type="evidence.txt"
    )


def log_transform_intensities(df: pd.DataFrame) -> pd.DataFrame:
    logger.info("Log10 transforming intensities")
    tmt_channels = utils.get_tmt_channels(df)
    df.loc[:, tmt_channels.columns] = np.log10(tmt_channels.replace(0, np.nan))
    return df


def sum_peptide_intensities(
    df: pd.DataFrame, debug: bool = False, run_lfq: bool = False
) -> pd.DataFrame:
    """sum up intensities per p-peptide across fractions and charge states"""
    logger.info(
        "Summing up intensities per p-peptide across fractions and charge states"
    )
    
    # Perform replacement
    df['Modified sequence'] = df["Modified sequence"].str.replace("M(Oxidation (M))", "M", regex=False)

    # TODO: split strings on semicolon before making unique
    def csv_list_unique(x):
        return ";".join(map(str, list(dict.fromkeys(x))))

    def aggregate_imputations(x):
        annotations = dict.fromkeys(x)
        if not set(annotations).issubset({"imputed;", ""}):
            raise ValueError(
                f"Found other annotations ({annotations}) besides imputations, need new solution to detect partially imputed peptides"
            )

        if set(annotations) == {"imputed;", ""}:
            return "partially imputed;"

        return ";".join(map(str, list(annotations)))

    # Phospho preprocessing
    # if 'Modified sequence' in df.columns and debug and run_lfq is False:
    #     print('hello')
    #     print(df.columns)
    #     df = df.groupby(['Batch', 'Modified sequence', 'Gene names'])
    #     print('step 1')
    #
    #     df = df.agg(
    #         # **{
    #         #     'new_genes': pd.NamedAgg(column='new_genes', aggfunc=csv_list_unique),
    #         # },
    #         **{f'Reporter intensity corrected {i}': pd.NamedAgg(column=f'Reporter intensity corrected {i}', aggfunc=sum) for i in
    #            range(1, 12)}
    #     ).reset_index()
    #     print('step 2')
    # Phospho raw for last module? TODO check this

    # TODO: handle the different cases more elegantly
    if "Modified sequence" in df.columns and run_lfq is False:
        df = df.groupby(["Batch", "Modified sequence"])

        # if 'Transferred spectra count' in df.columns:
        #     df = df.agg(
        #         **{
        #             'Proteins': pd.NamedAgg(column='Proteins', aggfunc=csv_list_unique),
        #             'Gene names': pd.NamedAgg(column='Gene names', aggfunc=csv_list_unique),
        #             'Transferred spectra count': pd.NamedAgg(column='Transferred spectra count', aggfunc=csv_list_unique),
        #         },
        #         **{f'Reporter intensity corrected {i}': pd.NamedAgg(column=f'Reporter intensity corrected {i}', aggfunc=sum) for i
        #            in
        #            range(1, 12)},
        #         **{f'Identification metadata {i}': pd.NamedAgg(column=f'Identification metadata {i}',
        #                                                        aggfunc=aggregate_imputations)
        #            for i in
        #            range(1, 12)}
        #     ).reset_index()
        # else:
        df = df.agg(
            **{
                "Proteins": pd.NamedAgg(column="Proteins", aggfunc=csv_list_unique),
                "Gene names": pd.NamedAgg(column="Gene names", aggfunc=csv_list_unique),
            },
            **{
                f"Reporter intensity corrected {i}": pd.NamedAgg(
                    column=f"Reporter intensity corrected {i}", aggfunc="sum"
                )
                for i in range(1, 12)
            },  # TODO: replace this hard coded number of TMT channels
            **{
                f"Identification metadata {i}": pd.NamedAgg(
                    column=f"Identification metadata {i}", aggfunc=aggregate_imputations
                )
                for i in range(1, 12)
            },  # TODO: replace this hard coded number of TMT channels
        ).reset_index()
    elif "Modified sequence" in df.columns and run_lfq:
        df = df.groupby(["Batch", "Modified sequence"])
        df = df.agg(
            **{
                "Proteins": pd.NamedAgg(column="Proteins", aggfunc=csv_list_unique),
                "Gene names": pd.NamedAgg(column="Gene names", aggfunc=csv_list_unique),
            },
            **{
                "Reporter intensity corrected 1": pd.NamedAgg(
                    column="Reporter intensity corrected 1", aggfunc="sum"
                )
            },
            **{
                "Identification metadata 1": pd.NamedAgg(
                    column="Identification metadata 1", aggfunc=aggregate_imputations
                )
            },
        ).reset_index()
    else:  # TODO: what is this case used for??
        df = df.groupby(["Batch", "Gene names"])
        df = df.agg(
            **{
                "Proteins": pd.NamedAgg(column="Proteins", aggfunc=csv_list_unique),
                "new_genes": pd.NamedAgg(column="new_genes", aggfunc=csv_list_unique),
            },
            **{
                f"Reporter intensity corrected {i}": pd.NamedAgg(
                    column=f"Reporter intensity corrected {i}", aggfunc="sum"
                )
                for i in range(1, 12)
            },  # TODO: replace this hard coded number of TMT channels
        ).reset_index()
    return df


def load_and_normalize(
    data_loader: "DataLoader",
    results_folder: Union[str, Path],
    sample_annotation_df: pd.DataFrame,
    data_type: str,
    normalize_to_reference: bool = False,
    debug: bool = True,
) -> pd.DataFrame:
    """
    Function for preprocessing MQ evidence files
    :param data_loader:
    :param results_folder: location of log and intermediately processed files
    :param data_type:
    :param debug: save extra output useful for debugging
    :return: Dataframe with combined preprocessed data from all batches
    """
    save_debug_df = get_save_debug_df_function(debug, results_folder, data_type)
    save_correction_factors = get_save_correction_factors_function(
        results_folder, data_type
    )

    # Parse evidence files and combine all batches
    all_batches = data_loader.load_data(MQ_EVIDENCE_COLUMNS)
    df = pd.concat(all_batches, ignore_index=True)
    save_debug_df(df, "_complete_raw")
    save_qc_info(df, sample_annotation_df, results_folder, data_type)

    # apply high dynamic range filter 
    logger.info(df.isna().sum().sum())
    df_new = df.copy()
    df = filter_high_dr_evidence(df, df_new)
    logger.info(df.isna().sum().sum())

    # MS2/MS3 intensity median centering within each batch
    df, correction_factors = data_loader.median_centering_within_batch(df)
    save_debug_df(df, "_after_1st_median")
    save_correction_factors(correction_factors, "_in_batch_correction_factors")

    # MS1 intensity median centering across batches
    df, correction_factors = data_loader.median_centering_ms1(df)
    save_debug_df(df, "_after_ms1_centering")
    save_correction_factors(correction_factors, "_ms1_correction_factors")

    ref_channels_df = sample_annotation_df.loc[
        sample_annotation_df["is_reference"] == True, :
    ]
    if normalize_to_reference:
        # Scale MS1 intensities such that reference channel MS1 contribution is constant across batches
        df = data_loader.scale_ms1_to_reference(df, ref_channels_df)
        save_debug_df(df, "_after_ms1_normalize_to_reference")
    else:
        # Impute an MS1 value for scans without MS1
        # This uses the assumption that the reference channel abundance should be constant across batches
        df = data_loader.impute_ms1_intensity(df, ref_channels_df)
        save_debug_df(df, "_after_ms1_imputation")

    # Distribute MS1 intensity over the TMT channels
    df = data_loader.redistribute_ms1_intensity(df)
    save_debug_df(df, "_after_ms1_correction")

    return df


def filter_high_dr_evidence(evidence_df, df_new, ratio_threshold: float = 0.01):
    """
    Filter out evidence with a log ratio below the specified threshold.
    
    Parameters:
    df (DataFrame): The input DataFrame containing evidence data.
    logratio (int): The log ratio threshold for filtering.
    
    Returns:
    DataFrame: Filtered DataFrame with high DR evidence.
    """
    pat_cols = [col for col in evidence_df.columns if col.startswith('Reporter intensity corrected')]

    # compute row-wise max
    row_max = evidence_df[pat_cols].max(axis=1)

    # compute ratio
    ratio_df = evidence_df[pat_cols].div(row_max, axis=0)

    # mask values smaller than threshold (1/100)
    small_mask = ratio_df < 0.01

    df_new[pat_cols] = evidence_df[pat_cols].mask(small_mask)
    
    return df_new


def get_save_debug_df_function(
    debug: bool, results_folder: str, data_type: str
) -> Callable[[pd.DataFrame, str], None]:
    def save_debug_df(df: pd.DataFrame, suffix: str) -> None:
        if not debug:
            return

        df.to_csv(
            os.path.join(results_folder, f"debug_preprocessed_{data_type}{suffix}.csv"),
            index=False,
            float_format="%.4g",
        )

    return save_debug_df


def get_save_correction_factors_function(
    results_folder: str, data_type: str
) -> Callable[[pd.DataFrame, str], None]:
    def save_correction_factors(df: pd.DataFrame, suffix: str) -> None:
        df.to_csv(
            os.path.join(results_folder, f"{data_type}{suffix}.csv"),
            index=True,
            float_format="%.4g",
        )

    return save_correction_factors


def save_qc_info(
    df: pd.DataFrame,
    sample_annotation_df: pd.DataFrame,
    results_folder: str,
    data_type: str,
):
    # Log QC info like identifications, summed intensity, median
    if data_type == "fp":
        retrieve_qc_info_function = retrieve_fp_qc_info
    elif data_type == "pp":
        retrieve_qc_info_function = retrieve_pp_qc_info
    else:
        raise ValueError(
            "Unknown data type for retrieve_qc_info: {data_type}. Valid values are 'fp' and 'pp'."
        )
    retrieve_qc_info_function(df, sample_annotation_df, results_folder, data_type)


def retrieve_fp_qc_info(df, sample_annotation_df, results_folder, data_type):
    print("Retrieving QC info for: ", data_type)

    df = filter_data(df, "fp")

    tmt_columns = df.loc[
        :, df.columns.str.startswith("Reporter intensity corrected")
    ].columns

    # sum up across fractions and charge states and subset to tmt columns
    intensities = get_tmt_intensities_per_peptide(df, tmt_columns)

    # Calculate count, sum, and mean for all columns at once using DataFrame methods
    peptide_count = intensities.groupby(["Batch"]).count().add_suffix("_Count")
    peptide_sum = intensities.groupby(["Batch"]).sum().add_suffix("_Sum")
    peptide_median = intensities.groupby(["Batch"]).median().add_suffix("_Median")

    result_df = pd.concat([peptide_count, peptide_sum, peptide_median], axis=1)

    melt_and_pivot_qc_info(result_df, sample_annotation_df, results_folder, data_type)
    get_batch_wise_qc_info(peptide_sum, peptide_median, results_folder, data_type)


def retrieve_pp_qc_info(df, sample_annotation_df, results_folder, data_type):
    print("Retrieving QC info for: ", data_type)
    df = filter_data(
        df, "fp"
    )  # filter only reverse and contaminants no matter data type
    # TODO: should this be "pp" instead of "fp"?

    tmt_columns = df.loc[
        :, df.columns.str.startswith("Reporter intensity corrected")
    ].columns

    phospho_only_df = df[df["Modifications"].apply(lambda x: "Phospho (STY)" in x)]
    phospho_only_df = phospho_only_df.groupby(["Batch", "Modified sequence"]).agg("sum")
    phospho_only_df.loc[:, tmt_columns] = phospho_only_df.loc[:, tmt_columns].replace(
        0.0, np.nan
    )

    # sum up across fractions and charge states and subset to tmt columns
    intensities = get_tmt_intensities_per_peptide(df, tmt_columns)

    # Calculate count, sum, and mean for all columns at once using DataFrame methods
    peptide_count = intensities.groupby(["Batch"]).count().add_suffix("_Count")

    intensities = phospho_only_df.filter(like="Reporter intensity corrected")
    pp_peptide_count = (
        intensities.groupby(["Batch"]).count().add_suffix("_Phospho-count")
    )
    peptide_sum = intensities.groupby(["Batch"]).sum().add_suffix("_Sum")
    peptide_median = intensities.groupby(["Batch"]).median().add_suffix("_Median")

    result_df = pd.concat(
        [peptide_count, pp_peptide_count, peptide_sum, peptide_median], axis=1
    )

    melt_and_pivot_qc_info(result_df, sample_annotation_df, results_folder, data_type)
    get_batch_wise_qc_info(peptide_sum, peptide_median, results_folder, data_type)


def get_tmt_intensities_per_peptide(df, tmt_columns):

    df = df.groupby(["Batch", "Modified sequence"]).agg("sum")
    df.loc[:, tmt_columns] = df.loc[:, tmt_columns].replace(0.0, np.nan)
    intensities = df.filter(like="Reporter intensity corrected")
    return intensities


def melt_and_pivot_qc_info(
    result_df: pd.DataFrame,
    sample_annotation_df: pd.DataFrame,
    results_folder: Union[str, Path],
    data_type: str,
):
    result_df.reset_index(inplace=True)
    melted_df = pd.melt(
        result_df, id_vars=["Batch"], var_name="Metric", value_name="Value"
    )
    melted_df[["Sample", "Metric_Type"]] = melted_df["Metric"].str.split(
        "_", expand=True
    )
    melted_df.drop(columns="Metric", inplace=True)
    final_df = melted_df.pivot_table(
        index=["Batch", "Sample"], columns="Metric_Type", values="Value"
    ).reset_index()

    final_df["Batch"] = final_df["Batch"].apply(lambda x: x.split("Batch")[1])
    final_df["Sample"] = final_df["Sample"].apply(lambda x: x.split(" ")[-1])

    sample_annotation_df["Batch Name"] = sample_annotation_df["Batch Name"].astype(str)
    sample_annotation_df["TMT Channel"] = sample_annotation_df["TMT Channel"].astype(
        str
    )

    merged_df = final_df.merge(
        sample_annotation_df[["Sample name", "Batch Name", "TMT Channel"]],
        left_on=["Batch", "Sample"],
        right_on=["Batch Name", "TMT Channel"],
        how="inner",
    )
    merged_df = merged_df[final_df.columns.tolist() + ["Sample name"]]
    merged_df.to_csv(
        os.path.join(results_folder, f"{data_type}_qc_numbers.csv"),
        index=True,
        float_format="%.4g",
    )


def get_batch_wise_qc_info(
    peptide_sum: pd.DataFrame,
    peptide_median: pd.DataFrame,
    results_folder: Union[str, Path],
    data_type: str,
):
    batch_sum = peptide_sum.sum(axis=1).add_suffix("_Sum")
    batch_median = peptide_median.median(axis=1).add_suffix("_Median")
    result_df = pd.concat([batch_sum, batch_median])
    result_df.to_csv(
        os.path.join(results_folder, f"{data_type}_qc_numbers_batch_wise.csv"),
        index=True,
        float_format="%.4g",
    )


def convert_long_to_wide_format(
    df: pd.DataFrame, has_metadata_cols: bool = False
) -> pd.DataFrame:
    """
    Converts dataframe from long format (11 reporter channels, batches concatenated row-wise)
    to wide format (samples as columns)
    """
    logger.info("Converting from long to wide format")

    tmt_channels = utils.get_tmt_channels(df).columns
    keep_columns = tmt_channels.tolist()

    if has_metadata_cols:
        metadata_cols = id_meta.as_metadata_columns(tmt_channels.str)
        keep_columns += metadata_cols.tolist()

    if "Modified sequence" in df.columns:
        keep_columns += ["Batch", "Gene names", "Proteins"]
        if "Transferred spectra count" in df.columns:
            keep_columns += ["Transferred spectra count"]

        index_col = "Modified sequence"
    else:
        keep_columns += ["Batch"]
        index_col = ["Gene names", "Proteins"]
    indexed_df = df.set_index(index_col)[keep_columns]

    all_batches = list()
    for batch_name, df in indexed_df.groupby("Batch"):
        df = df.drop(columns="Batch")
        df = df.rename(columns=lambda x: f"{x} {batch_name}")
        all_batches.append(df)
    wide_df = pd.DataFrame().join(all_batches, how="outer")
    wide_df = wide_df.reset_index()

    if has_metadata_cols:
        wide_df = aggregate_csv_columns(wide_df, "Gene names")
        wide_df = aggregate_csv_columns(wide_df, "Proteins")
    return wide_df


def aggregate_csv_columns(df: pd.DataFrame, column_name: str) -> pd.DataFrame:
    all_cols = df.columns[df.columns.str.startswith(column_name)]

    def csv_list_unique(csv_string):
        l = csv_string.split(";")
        s = sorted({x for x in l if x != "nan"})
        return ";".join(s)

    df[column_name] = (
        df[all_cols].astype(str).agg(";".join, axis=1).apply(csv_list_unique)
    )

    df = df.drop(columns=all_cols)
    return df


def drop_duplicate_indices(df: pd.DataFrame) -> pd.DataFrame:
    """Keeps only first entry for entries with duplicated indices"""
    return df[~df.index.duplicated(keep="first")]


def merge_ms1_columns(df: pd.DataFrame) -> pd.DataFrame:
    ms1s = list()
    for _, df in df.groupby("Batch"):
        batch_name = df["Batch"].iloc[0]
        ms1_intensity_df = drop_duplicate_indices(
            df.set_index("Modified sequence")[["Intensity"]].rename(
                columns={"Intensity": batch_name}
            )
        )
        ms1s.append(ms1_intensity_df)

    # join MS1 intensities by modified sequence
    return pd.DataFrame().join(ms1s, how="outer")


def median_centering(df: pd.DataFrame) -> pd.DataFrame:
    """
    Normalizes samples by multiplying each sample with its own correction factor
    calculated by division of average median of reference channels with the median of
    the sample. That way sample distributions gets centered around the same value
    :param df: dataframe of batch intensities (11 reporter channels) or dict of MS1 values per batch
    """
    # Multiplying with correction factor to center samples
    df = df.replace(0, np.nan)
    medians = df[df.isna().sum(axis=1) <= 3].median(axis=0)
    correction_factor = medians.mean() / medians
    df = df * correction_factor
    return df, correction_factor


def impute_data(df: pd.DataFrame) -> pd.DataFrame:
    """
    If a peptide is not detected for a channel but is detected in the batch,
    impute with 1/100th of the maximum intensity or the minimum intensity for
    a peptide within its batch, whichever is lower.

    :param df: dataframe of batch intensities (11 reporter channels)
    """
    logger.info("Imputing data")

    tmt_cols_df = utils.get_tmt_channels(df)
    tmt_cols_df = tmt_cols_df.replace(0, np.nan)

    patient_channels_df = tmt_cols_df.filter(
        regex=r"^Reporter intensity corrected [1-8]$"
    )

    x = tmt_cols_df.max(axis=1) / 100
    y = tmt_cols_df.min(axis=1)
    imputation_values = pd.concat([x, y], axis="columns").min(axis="columns")

    # apply per column
    patient_channels = patient_channels_df.columns
    metadata_cols = id_meta.as_metadata_columns(patient_channels.str)

    def add_imputation_status(z):
        return np.where(z.isna() & imputation_values.notna(), "imputed;", "")

    df[metadata_cols] += tmt_cols_df.apply(add_imputation_status).rename(
        columns=id_meta.as_metadata_columns
    )
    df.loc[:, patient_channels] = patient_channels_df.apply(
        lambda z: z.fillna(imputation_values)
    )
    return df


def filter_data(df: pd.DataFrame, data_type: str) -> pd.DataFrame:
    """
    filters out contaminants, reverse sequences from decoy database,
    For phospho, non-phospho modified peptides are removed as well.
    :param df:
    :param results_folder:
    :param data_type:
    """
    logger.info(f"Filtering data {data_type}")

    logger.info(f"- before filtering: {df.shape[0]}")

    df = df[df["Potential contaminant"] != "+"]
    logger.info(f"- after contaminant removal: {df.shape[0]}")

    df = df[df["Reverse"] != "+"]
    logger.info(f"- after reverse removal: {df.shape[0]}")

    df = df.drop(["Reverse", "Potential contaminant"], axis=1)

    if data_type == "pp":
        df = df[df["Modifications"].str.contains("Phospho (STY)", regex=False)]
        logger.info(f"- after unmodified peptides removal: {df.shape[0]}")

        df = df.drop("Modifications", axis=1)

    return df


def rename_columns_with_sample_ids(
    df: pd.DataFrame,
    channel_to_sample_id_dict: Dict[str, str],
    index_cols: List[str],
) -> pd.DataFrame:
    """
    Transform column names of the format 'Reporter intensity corrected <TMT_channel> <batch_name>' to the sample names
    """
    df = df.set_index(index_cols)
    tmt_channels = list(channel_to_sample_id_dict.keys())
    metadata_cols = list(map(id_meta.as_metadata_columns, tmt_channels))

    keep_cols = tmt_channels + metadata_cols

    # if not remove_ref:
    #     ref_channels_cols = df.filter(regex='^Reporter intensity corrected (9|10|11)').columns.tolist()
    #     keep_cols += ref_channels_cols
    df = df.filter(items=keep_cols, axis=1)  # here it happens

    # build dictionary to also rename the metadata columns with the sample ids
    metadata_cols_with_sample_id = map(
        lambda x: f"Identification metadata {x}", channel_to_sample_id_dict.values()
    )
    metadata_to_sample_id_dict = dict(zip(metadata_cols, metadata_cols_with_sample_id))

    # add "pat_" prefix to patient intensity columns
    channel_to_sample_id_dict_final = {
        k: "pat_" + v
        for k, v in channel_to_sample_id_dict.items()
        if not v.startswith("ref_")
    }
    # keep "ref_" prefix for reference intensity columns
    channel_to_sample_id_dict_final |= {
        k: v for k, v in channel_to_sample_id_dict.items() if v.startswith("ref_")
    }

    rename_dict = {**channel_to_sample_id_dict_final, **metadata_to_sample_id_dict}
    df = df.rename(columns=rename_dict)
    df = df.reset_index()
    return df


def get_data_location(
    maxquant_super_folder: Union[str, Path], data_type: str, file_type="evidence.txt"
) -> Tuple[List[str], List[str]]:
    """
    Return file paths of all relevant files from maxquant_super_folder folder
    :param maxquant_super_folder: Folder containing results from MQ searches
    :param data_type:
    :param file_type:
    :return: List of file paths for full and phospho proteome data
    """
    if not os.path.isdir(maxquant_super_folder):
        raise IOError(
            f"MaxQuant super folder is not a directory: {maxquant_super_folder}"
        )

    validity_check = is_valid_fp_file
    if data_type == "pp":
        validity_check = is_valid_pp_file

    evidence_files = []
    for file_to_check in sorted(glob(
        os.path.join(maxquant_super_folder, "**", file_type), recursive=True
    )):
        if validity_check(file_to_check):
            evidence_files.append(file_to_check)
    evidence_files = sorted(evidence_files)
    return evidence_files


def filter_evidence_files(
    files: List[Union[str, Path]], data_type: str, batches: List[str]
):
    # TODO: test that we can crash it if multiple batch 40
    return [
        file
        for file in files
        for batch_number, cohort in batches
        if os.path.join(cohort, f"Batch{batch_number}_{data_type}") in file
    ]


def is_valid_fp_file(directory: Union[str, Path]) -> bool:
    return (
        "FP" in directory
        and "Merged_FP" not in directory
        and "QC_FAILED" not in directory
        and "outdated" not in directory
        and "Technician" not in directory
        and "_var_" not in directory
    )


def is_valid_pp_file(directory: Union[str, Path]) -> bool:
    return (
        "PP" in directory
        and "QC_FAILED" not in directory
        and "outdated" not in directory
        and "Technician" not in directory
        and "_var_" not in directory
    )
