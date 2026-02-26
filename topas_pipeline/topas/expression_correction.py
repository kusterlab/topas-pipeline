import collections
import sys
from pathlib import Path
import logging

import pandas as pd
import numpy as np
import scipy

from tqdm import tqdm
from ..preprocess import phospho_grouping
from .. import utils

# hacky way to get the package logger instead of just __main__ when running as a module
logger = logging.getLogger(__package__ + "." + Path(__file__).stem)


def correct_phospho_for_protein_expression(
    results_folder: str, sample_annotation_file: str, overwrite: bool = False
):
    results_folder = Path(results_folder)
    expression_corrected_file = (
        results_folder / "preprocessed_pp2_agg_batchcorrected_expressioncorrected.csv"
    )
    if expression_corrected_file.is_file():
        if not overwrite:
            logger.info(
                f"Reusing previously generated expression corrected intensities: {expression_corrected_file}"
            )
            return
        logger.info(f"Found existing results but overwrite flag was set.")

    phospho_df = phospho_grouping.read_cohort_intensities_df(
        results_folder / "preprocessed_pp2_agg_batchcorrected.csv",
        sample_annotation_file=sample_annotation_file,
        keep_identification_metadata_columns=False,
    )
    full_df = load_full_proteome_df(results_folder)

    # only use patient channels for predictions
    patient_columns = utils.filter_for_patient_columns(phospho_df).columns

    if not all(phospho_df.columns.sort_values() == full_df.columns.sort_values()):
        raise ValueError(
            f"Missing columns in full: {set(phospho_df.columns) - set(full_df.columns)}\nMissing columns in phospho: {set(full_df.columns) - set(phospho_df.columns)}"
        )

    mod_seq_df = get_corresponding_protein_expression(
        phospho_df, full_df, patient_columns
    )

    # predict full proteome from phospho to fill in missing protein expressions
    full_predicted = predict_full_from_phospho(phospho_df, mod_seq_df, patient_columns)

    # correct phospho for protein expression
    phospho_corrected = predict_phospho_from_full(
        phospho_df,
        mod_seq_df,
        full_predicted,
        patient_columns,
    )

    logger.info(
        f"Writing expression corrected phospho matrix to {expression_corrected_file}"
    )
    phospho_corrected.to_csv(expression_corrected_file, float_format="%.4f")


def load_full_proteome_df(results_folder: str) -> pd.DataFrame:
    logger.info("Loading full proteome data")

    header = pd.read_csv(results_folder / "preprocessed_fp.csv", index_col=0, nrows=1)
    intensity_cols = utils.filter_for_sample_columns(header).columns.tolist()

    dtype_dict = collections.defaultdict(lambda: "str")
    dtype_dict |= {c: "float32" for c in intensity_cols}
    full_df = pd.read_csv(
        results_folder / "preprocessed_fp.csv",
        index_col=0,
        usecols=intensity_cols + ["Gene names"],
        dtype=dtype_dict,
    )
    full_df.index = full_df.index.fillna("")
    return full_df


def get_corresponding_protein_expression(
    phospho_df: pd.DataFrame,
    full_df: pd.DataFrame,
    patient_columns: pd.Index,
) -> pd.DataFrame:
    # sum multiple gene expressions if the same phospho peptide matches to multiple genes

    # expand for single gene names matching
    mod_seq_df = phospho_df.index.to_frame().reset_index(drop=True)
    mod_seq_df["Gene names"] = mod_seq_df["Gene names"].str.split(";")
    mod_seq_df = mod_seq_df.explode(column="Gene names")

    # back aggregate to unique phospho
    mod_seq_df = mod_seq_df.merge(right=10**full_df, on="Gene names", how="left")
    mod_seq_df = mod_seq_df.merge(
        right=phospho_df.index.to_frame().reset_index(drop=True),
        on=phospho_df.index.names,
        how="left",
    )

    # divide into numeric, string and groupby col and create agg dict
    numeric_cols = mod_seq_df.select_dtypes(include=np.number).columns
    group_col = "Modified sequence group"
    string_cols = [
        c for c in phospho_df.index.names if c != group_col and c not in numeric_cols
    ]
    agg_dict = {
        col: lambda x: ";".join(x.dropna().astype(str).unique()) for col in string_cols
    }
    for col in numeric_cols:
        agg_dict[col] = "sum"

    # Group by and aggregate --> then log10 transform
    aggregated_df = mod_seq_df.groupby("Modified sequence group", as_index=False).agg(
        agg_dict
    )
    aggregated_df[numeric_cols] = aggregated_df[numeric_cols].replace(0, np.nan)
    aggregated_df[numeric_cols] = np.log10(aggregated_df[numeric_cols])
    aggregated_df = aggregated_df.set_index(phospho_df.index.names)

    # mask low observations to prevent false corrections
    aggregated_df[(aggregated_df[patient_columns].isna().mean(axis=1) > 0.8)] = 0

    return aggregated_df


def predict_full_from_phospho(
    phospho_df: pd.DataFrame, mod_seq_df: pd.DataFrame, patient_columns: pd.Index
) -> pd.DataFrame:
    logger.info("Predicting protein expression from phospho proteome")
    # for fast analysis switch to numpy
    P = phospho_df[patient_columns].values
    F = mod_seq_df.loc[phospho_df.index, patient_columns].values

    # Predict the fullproteome based on phospho
    linear_fits = {}
    for i, idx in tqdm(enumerate(phospho_df.index)):
        linear_fits[idx] = [*get_bounded_model(x=P[i], y=F[i])]
    linear_fits = pd.DataFrame(
        linear_fits, index=["N data points", "slope", "intersept", "var ratio"]
    ).T
    linear_fits.index = linear_fits.index.set_names(phospho_df.index.names)

    # Predict full based on phospho
    slopes = linear_fits.slope.values
    intercepts = linear_fits.intersept.values
    full_predicted = pd.DataFrame(
        (phospho_df.T.values * slopes + intercepts).T,
        columns=phospho_df.columns,
        index=linear_fits.index,
    )
    return full_predicted


def predict_phospho_from_full(
    phospho_df: pd.DataFrame,
    mod_seq_df: pd.DataFrame,
    full_predicted: pd.DataFrame,
    patient_columns: pd.Index,
) -> pd.DataFrame:
    logger.info("Predicting phospho from protein expression")
    # for fast analysis switch to numpy
    P = phospho_df[patient_columns].values
    F = mod_seq_df.loc[phospho_df.index, patient_columns].values

    # Predict the phospho based on fullproteome -> IF phospho is predictable (this amount needs to be subtracted)
    linear_fits = {}
    for i, idx in tqdm(enumerate(phospho_df.index)):
        linear_fits[idx] = [*get_bounded_model(x=F[i], y=P[i])]

    linear_fits = pd.DataFrame(
        linear_fits, index=["N data points", "slope", "intersept", "var ratio"]
    ).T
    linear_fits.index = linear_fits.index.set_names(phospho_df.index.names)

    # Impute fullprotome NaNs with predicted values
    F_imputed = mod_seq_df.loc[phospho_df.index, phospho_df.columns].values
    F_pred = full_predicted.loc[phospho_df.index, phospho_df.columns].values
    F_imputed[np.isnan(F_imputed)] = F_pred[np.isnan(F_imputed)]

    # Clip low abundant outlier to 2.5% quantile to prevent over correction
    F_q25 = mod_seq_df[patient_columns].quantile(q=0.025, axis=1).values
    F_imputed = np.clip(F_imputed.T, a_min=F_q25, a_max=None).T
    F = np.clip(F.T, a_min=F_q25, a_max=None).T

    # Correct values
    slopes = linear_fits.slope.replace(np.nan, 0).values
    intercepts = linear_fits.intersept.replace(np.nan, 0).values

    # corrected values = raw values - (fullproteome contributions) + (pushback to avg hights)
    phospho_corrected = (
        phospho_df.T.values
        - (F_imputed.T * slopes + intercepts)
        + (np.nanmean(F.T, axis=0) * slopes + intercepts)
    ).T
    phospho_corrected = pd.DataFrame(
        phospho_corrected, columns=phospho_df.columns, index=linear_fits.index
    )

    return phospho_corrected


def get_std_ratio(x, y, limits=(0.1, 0.1)):
    x_tstd = scipy.stats.mstats.trimmed_std(x, limits=limits)
    y_tstd = scipy.stats.mstats.trimmed_std(y, limits=limits)
    return y_tstd / x_tstd


def get_spread(x, limits=(0.0, 0.0)):
    spread = scipy.stats.mstats.trimmed_std(x, limits=limits) * 4
    return spread


def calculate_slope(x, y):
    x = x - x.mean()
    y = y - y.mean()
    return np.sum(x * y) / np.sum(x**2)


def calculate_intercept(x, y, slope):
    return y.mean() - slope * x.mean()


def get_bounded_model(x, y):
    mask = np.isfinite(x) & np.isfinite(y)
    n = mask.sum()

    if n:
        x, y = x[mask], y[mask]
        ratio = get_std_ratio(x, y)
        x_spread = get_spread(x)
        slope = calculate_slope(x, y)
        slope = np.clip(slope, 0, min(x_spread, ratio, 1))
        intercept = calculate_intercept(x, y, slope)
        return n, slope, intercept, ratio
    return 0, 0.0, 0.0, np.nan


"""
python3 -m topas_pipeline.topas.expression_correction -c config_patients.json
"""
if __name__ == "__main__":
    import argparse

    from .. import config

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-c", "--config", required=True, help="Absolute path to configuration file."
    )
    parser.add_argument(
        "-o",
        "--overwrite",
        action="store_true",
        help="Ignore existing results and recompute outputs.",
    )
    args = parser.parse_args(sys.argv[1:])

    configs = config.load(args.config)

    correct_phospho_for_protein_expression(
        results_folder=configs.results_folder,
        sample_annotation_file=configs.sample_annotation,
        overwrite=args.overwrite,
    )
