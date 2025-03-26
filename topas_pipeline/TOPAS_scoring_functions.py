import os
import random
import logging

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.stats import norm


logger = logging.getLogger(__name__)


def calculate_psite_weights(df: pd.DataFrame):
    """Calculates weight for each p-site according to the peptide occurrence.

    This function deals with the situation that a single p-site can be detected in
    multiple modified peptides. In the end, we want each p-site to only be weighted
    once, to prevent a single p-site from disproportionally contributing to the final
    TOPAS score.

    A peptide's occurrence is the number of patients the peptide was detected in. We
    trust the peptide (containing the p-site) with a higher occurrence more than those
    which were only detected sporadically.

    Args:
        df: DataFrame with z-scores in columns named 'pat_<patient_id>'

    Returns:
        Input dataframe with weight columns for each patient named 'weight_<patient_id>'
    """
    df = df.sort_values("Modified sequence")

    patient_zscore_columns = [col for col in df.columns if "pat_" in col]
    patient_weight_columns = [
        col.replace("pat_", "weight_") for col in patient_zscore_columns
    ]

    # calculate peptide count per patient, setting it to 0 where patient intensity is 0
    weights_df = df[patient_zscore_columns + ["Site positions"]]
    weights_df = weights_df.rename(columns=lambda x: x.replace("pat_", "weight_"))
    weights_df[patient_weight_columns] = (
        weights_df[patient_weight_columns].notna().mul(df["Peptide count"], axis=0)
    )

    # calculate psite count per patient
    psite_counts_df = weights_df.groupby("Site positions").sum()

    # divide peptide count by psite count per patient
    # TODO: check w unittest if the two steps also work
    # weights_df[patient_weight_columns] /= weights_df[['Site positions']].merge(psite_counts_df, on='Site positions', how='left', validate="many_to_one")[patient_weight_columns].values

    merged_df = weights_df[["Site positions"]].merge(
        psite_counts_df, on="Site positions", how="left", validate="many_to_one"
    )
    weights_df[patient_weight_columns] /= merged_df[patient_weight_columns].values

    weights_df = weights_df.drop(columns=["Site positions"])
    weights_df = weights_df.replace(0, np.nan)

    return pd.concat([df, weights_df], axis=1)


def calculate_modified_sequence_weights(
    patient_dataframe: pd.DataFrame, summing_column: str
):
    """Sums the p-site weights of each modified sequence constituting p-sites.

    Args:
        patient_dataframe: dataframe with p-site weights and z-scores
        summing_column: column name with the annotation that will eventually be grouped by

    Returns:
        dataframe with weights for each modified sequence
    """
    weight_dataframe = patient_dataframe.groupby([summing_column, "Modified sequence"])
    weight_dataframe = weight_dataframe.agg(
        **{
            weight_col: pd.NamedAgg(column=weight_col, aggfunc="sum")
            for weight_col in patient_dataframe.columns.tolist()
            if "weight_" in weight_col
        },
        **{
            pat_col: pd.NamedAgg(column=pat_col, aggfunc="first")
            for pat_col in patient_dataframe.columns.tolist()
            if "pat_" in pat_col
        },
    ).reset_index()
    return weight_dataframe


def cap_zscores_and_weights(patient_dataframe):
    """Cap all z-scores between -4 and 4 and all weights between 0 and 1.

    Args:
        patient_dataframe: dataframe with modified sequence weights and z-scores

    Returns:
        dataframe with capped modified sequence weights and z-scores
    """
    # Reuse for CJ pipeline FP and PP clipping?
    capped_dataframe = patient_dataframe.copy()
    weightcols = [col for col in capped_dataframe.columns if "weight_" in col]
    patcols = [col for col in capped_dataframe.columns if "pat_" in col]
    capped_dataframe[weightcols] = (
        capped_dataframe[weightcols].astype(float).clip(upper=1)
    )
    # capped_dataframe[weightcols] = capped_dataframe[weightcols].astype(float)
    # TODO: Test if this capping step is required too
    capped_dataframe[patcols] = (
        capped_dataframe[patcols].astype(float).clip(upper=4, lower=-4)
    )
    return capped_dataframe


def calculate_weighted_z_scores(df: pd.DataFrame):
    """Multiplies the patient z-scores by the p-site weight

    Args:
        df: DataFrame with z-scores in columns named 'pat_<patient_id>' and weights in columns named 'weight_<patient_id>'

    Returns:
        Input dataframe with added columns for each patient named 'weighted_<patient_id>'
    """
    patient_zscore_columns = [col for col in df.columns if "pat_" in col]
    patient_weight_columns = [
        col.replace("pat_", "weight_") for col in patient_zscore_columns
    ]

    weighted_zscores_df = df[patient_zscore_columns].mul(
        df[patient_weight_columns].rename(
            columns=lambda x: x.replace("weight_", "pat_")
        )
    )
    weighted_zscores_df = weighted_zscores_df.rename(
        columns=lambda x: x.replace("pat_", "weighted_")
    )

    return pd.concat([df, weighted_zscores_df], axis=1)


def sum_weighted_z_scores(patient_dataframe, by):
    relevant_columns = [
        col
        for col in patient_dataframe.columns
        if ("weight_" not in col) and ("pat_" not in col)
    ]
    score_dataframe = patient_dataframe[relevant_columns].groupby([by])
    # TODO: move aggregation into var for readability
    score_dataframe = score_dataframe.agg(
        **(
            {
                score: pd.NamedAgg(column=score, aggfunc="sum")
                for score in patient_dataframe.columns
                if "weighted_" in score
            }
        )
    ).reset_index()
    score_dataframe = score_dataframe.replace(0, np.nan)
    score_dataframe = score_dataframe.rename(
        columns={
            col: col.replace("weighted_", "pat_") for col in score_dataframe.columns
        }
    )
    return score_dataframe


def plot_histograms_to_check_normality(score_dataframe, number_of_plots=20):

    for row in random.sample(range(len(score_dataframe)), number_of_plots):
        series = score_dataframe.iloc[row]
        zscores = series[1:-2].astype(float)
        if zscores.isna().all():
            logger.info(f"{series[0]} only contains NaNs! Not plotting...")
            continue
        mean, std = norm.fit(zscores)
        zscores.plot.density(label="Z-score PDF")
        zscores = np.array(zscores)
        plt.hist(zscores, weights=np.zeros_like(zscores) + 1.0 / zscores.size)
        xmin, xmax = plt.xlim()
        x = np.linspace(xmin, xmax, 100)
        p = norm.pdf(x, mean, std)
        plt.plot(x, p, "r", linewidth=2, label="Normal distribution")
        plt.axvline(x=mean, linestyle="-", color="r")
        plt.legend()
        plt.title(series[0])
        plt.show()
        plt.close()


def second_level_z_scoring(patient_dataframe, by_column, plot_histograms=False):
    # TODO: refactor this code
    patient_dataframe["mean"] = patient_dataframe[
        [
            col
            for col in patient_dataframe.columns
            if col not in [by_column, "Modified sequence"]
        ]
    ].mean(axis=1)
    patient_dataframe["stdev"] = patient_dataframe[
        [
            col
            for col in patient_dataframe.columns
            if col not in [by_column, "Modified sequence", "mean"]
        ]
    ].std(axis=1)
    for patient in patient_dataframe.columns:
        if patient in [by_column, "Modified sequence", "mean", "stdev"]:
            continue
        patient_dataframe[patient] = (
            patient_dataframe[patient] - patient_dataframe["mean"]
        ) / patient_dataframe["stdev"]
    if plot_histograms:
        plot_histograms_to_check_normality(patient_dataframe)
    return patient_dataframe


def calculate_peptide_occurrence(pp_df: pd.DataFrame):
    pp_df = pp_df.replace("nan", np.nan)
    pp_df = pp_df.replace("", np.nan)
    patient_list = pp_df.filter(regex="pat_").columns.tolist()
    pp_df = pp_df.astype({patient: "float32" for patient in patient_list})
    # pp_df = pp_df.where(pp_df[patient_list] < 4, 4)
    # pp_df = pp_df.where(pp_df[patient_list] > -4, -4)
    pp_df["Peptide count"] = pp_df[patient_list].notna().sum(axis=1)
    pp_df["Peptide occurrence"] = (
        pp_df["Peptide count"].astype(str) + "/" + str(len(patient_list))
    )
    return pp_df


def topas_score_preprocess(results_folder, discard_isoforms=True):

    # TODO: make a wrapper/decorator that checks if exists and if not run this function (use when this function is run from other modules)

    filepath = os.path.join(results_folder, "topas_score_preprocessed.tsv")

    if os.path.exists(filepath):
        logger.info(f"Reading previously generated file: {filepath}")
        return pd.read_csv(filepath, sep="\t")

    patients_proteins = pd.read_csv(
        os.path.join(results_folder, "annot_pp.csv"), index_col=0
    )
    patients_proteins = patients_proteins.drop(
        patients_proteins.filter(regex="pat_|ref_").columns, axis=1
    )

    patients_zscores = pd.read_csv(
        os.path.join(results_folder, "phospho_measures_z.tsv"),
        keep_default_na=False,
        sep="\t",
    )
    patients_zscores = patients_zscores.rename(
        columns={
            colname: colname.replace("zscore_", "")
            for colname in patients_zscores.columns
        }
    )

    patient_columns = patients_zscores.filter(regex="pat_").columns.tolist()
    patients_zscores = patients_zscores.loc[
        :, ["Gene names", "Modified sequence"] + patient_columns
    ]

    drop_cols = [
        "PSP_URL",
        "PSP_URL_extra",
        "basket",
        "basket_weights",
        "sub_basket",
        "sub_basket_weights",
        "other",
        "other_weights",
        "rtk",
        "rtk_weights",
        "drug",
        "drug_weights",
        "Proteins",
    ]
    patients_proteins = patients_proteins.drop(columns=drop_cols, errors="ignore")
    metadata_cols = patients_proteins.filter(
        regex="^Identification metadata "
    ).columns.tolist()
    patients_proteins = patients_proteins.drop(columns=metadata_cols, errors="ignore")

    patients = pd.merge(
        left=patients_zscores,
        right=patients_proteins,
        on=["Modified sequence", "Gene names"],
        how="inner",
        validate="one_to_one",
    )
    patients.set_index("Modified sequence", inplace=True)

    # for debugging: if using preprocessed_pp instead of annot_pp, annotate the p-sites here
    # clinical_tools.add_phospho_annotations(patients, pspFastaFile, pspKinaseSubstrateFile, pspAnnotationFile, pspRegulatoryFile)

    patients = calculate_peptide_occurrence(patients)  # better name?

    patients.rename(
        columns={
            "Site positions": "All site positions",
            "Site positions identified (MQ)": "Site positions",
        },
        inplace=True,
    )

    patients = patients[patients["Site positions"].str.len() > 0]

    patients["Site positions"] = patients["Site positions"].str.split(";")
    patients = patients.explode("Site positions")
    if discard_isoforms:
        # TODO: comment with regex example / documentation or use function for this
        patients = patients[
            ~patients["Site positions"].str.contains(r"^(?:[^\W_]+-\d+_[STY]\d+)$")
        ]

    patients = patients.reset_index()
    patients.to_csv(filepath, sep="\t", float_format="%.4g")
    return patients


def read_preprocessed_df(location: str):
    df = pd.read_csv(location, sep="\t")
    # maybe set index col?
    return df


def calculate_per_patient_targets(scored_peptide_df: pd.DataFrame, grouping_by: str):
    patients = [i for i in scored_peptide_df.columns if "pat_" in i]
    per_patient_targets = (
        scored_peptide_df.groupby(grouping_by)[patients].count().reset_index()
    )
    per_patient_targets.rename(
        columns={
            col: col.replace("pat_", "targets_")
            for col in per_patient_targets.columns
            if "pat_" in col
        },
        inplace=True,
    )
    return per_patient_targets


def get_target_space(annotated_peptides_df, grouping_by, scored_peptides_df):
    total_targets = (
        annotated_peptides_df.groupby(grouping_by)["Modified sequence"]
        .nunique()
        .reset_index(name="No. of total targets")
    )
    per_patient_targets = calculate_per_patient_targets(scored_peptides_df, grouping_by)
    target_space = pd.merge(
        left=total_targets,
        right=per_patient_targets,
        on=grouping_by,
        how="inner",
        validate="1:1",
    )
    return target_space


if __name__ == "__main__":
    print("This is not executable!")
