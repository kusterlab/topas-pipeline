import sys
import logging
from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# hacky way to get the package logger instead of just __main__ when running as module
logger = logging.getLogger(__package__ + "." + __file__)


PLOT_Y_LIMITS = {
    "FP": {"MS scans": 2500, "MS/MS scans": 10000, "ID rate": 70, "PSMs": 7000},
    "PP": {"MS scans": 6000, "MS/MS scans": 20000, "ID rate": 70, "PSMs": 7000},
}


def main(mq_folder_base_path: Path, search_type: str):
    if "combined" in str(mq_folder_base_path):
        mq_folder_base_path = mq_folder_base_path.parents[2]
    mq_folder_path = mq_folder_base_path / "combined" / "txt"

    if search_type not in ["FP", "PP"]:
        raise ValueError(f"Search type {search_type} is not one of FP or PP.")

    # write results to subfolder of original txt folder
    qc_plot_output_path = mq_folder_path / "QC_plots"

    # skip if folder already contains output table
    if (qc_plot_output_path / f"ID_table_{search_type}.csv").is_file():
        # logger.info(
        #     f"Skipped generating QC plots, found existing results in {qc_plot_output_path}"
        # )
        return

    logger.info(f"Generating QC plots for {mq_folder_base_path.name}")
    qc_plot_output_path.mkdir(exist_ok=True)

    # use relative paths to find tracking table and basket table paths
    tracking_table_directory_path = mq_folder_base_path.parent

    # QC plots for summary.txt
    summary_df, summary_df_totals = read_summary_txt(mq_folder_path)
    plot_summary_qc_plots(
        summary_df, summary_df_totals, qc_plot_output_path, search_type
    )

    # QC plots for evidence.txt
    evidence_df = read_evidence_txt(mq_folder_path)
    batch_name = evidence_df.iloc[0]["Experiment"]
    peptide_intensities_per_fraction, peptide_intensities_per_channel = (
        plot_evidence_qc_plots(
            evidence_df, qc_plot_output_path, batch_name, prefix="Peptides"
        )
    )

    ppeptide_intensities_per_fraction = None
    ppeptide_intensities_per_channel = None
    imac_selectivity_df = None
    if search_type == "PP":
        evidence_df_only_phosho = evidence_df[evidence_df["Phospho (STY)"] > 0]
        ppeptide_intensities_per_fraction, ppeptide_intensities_per_channel = (
            plot_evidence_qc_plots(
                evidence_df_only_phosho,
                qc_plot_output_path,
                batch_name,
                prefix="Ppeptides",
            )
        )

        imac_selectivity_df = plot_imac_selectivity_per_fraction(
            evidence_df, evidence_df_only_phosho, qc_plot_output_path, batch_name
        )

    # QC plots for proteinGroups.txt
    if search_type == "FP":
        protein_groups_df = read_protein_groups_txt(mq_folder_path)
    elif search_type == "PP":
        protein_groups_df = group_evidence_by_gene_name(evidence_df)

    basket_table_path = tracking_table_directory_path.parent / "baskets.txt"
    if basket_table_path.is_file():
        basket_members = read_basket_annotation_file(basket_table_path)
        plot_basket_identification_rates(
            basket_members, protein_groups_df, qc_plot_output_path
        )

    # Table with per-fraction statistics
    write_batch_statistics_table(
        batch_name,
        peptide_intensities_per_fraction,
        peptide_intensities_per_channel,
        ppeptide_intensities_per_fraction,
        ppeptide_intensities_per_channel,
        imac_selectivity_df,
        protein_groups_df,
        summary_df_totals,
        summary_df,
        qc_plot_output_path,
        search_type,
    )

    # Append per-channel statistics to tracking table across all batches and create plots as pdf
    tracking_table_intensities = peptide_intensities_per_channel
    if search_type == "PP":
        tracking_table_intensities = ppeptide_intensities_per_channel

    tracking_table = append_channel_statistics_to_tracking_table(
        batch_name,
        tracking_table_intensities,
        protein_groups_df,
        tracking_table_directory_path,
        search_type,
    )
    plot_tracking_table_pdf(tracking_table, tracking_table_directory_path, search_type)


def read_summary_txt(mq_folder_path: Path):
    summary_df = pd.read_csv(mq_folder_path / "summary.txt", sep="\t")
    summary_df = summary_df.rename(
        columns={
            "MS": "MS scans",
            "MS/MS": "MS/MS scans",
            "MS/MS Identified [%]": "ID rate",
            "Peptide Sequences Identified": "PSMs",
        }
    )
    # for new MQ version
    summary_df = summary_df.rename(
        columns={
            "MS": "MS scans",
            "MS/MS": "MS/MS scans",
            "MS/MS identified [%]": "ID rate",
            "Peptide sequences identified": "PSMs",
        }
    )

    # Store row with column totals calculated by MQ in separate dataframe
    summary_df_totals = summary_df.iloc[-1]

    # Store summary dataframe without row with column totals
    summary_df = summary_df.iloc[:-1]
    summary_df["Fraction"] = summary_df["Fraction"].astype(int)

    return summary_df, summary_df_totals


def plot_summary_qc_plots(
    summary_df: pd.DataFrame,
    summary_df_totals: pd.Series,
    output_path: Path,
    search_type: str,
):

    plot_summary_by_fraction(
        summary_df,
        output_path,
        column="MS scans",
        ylim=PLOT_Y_LIMITS[search_type]["MS scans"],
    )
    plot_summary_by_fraction(
        summary_df,
        output_path,
        column="MS/MS scans",
        ylim=PLOT_Y_LIMITS[search_type]["MS/MS scans"],
    )
    plot_summary_by_fraction(
        summary_df,
        output_path,
        column="ID rate",
        ylim=PLOT_Y_LIMITS[search_type]["ID rate"],
    )
    plot_summary_by_fraction(
        summary_df,
        output_path,
        column="PSMs",
        ylim=PLOT_Y_LIMITS[search_type]["PSMs"],
        title_suffix=f"\nTotal {summary_df_totals['PSMs']}",
    )


def read_evidence_txt(mq_folder_path: Path):
    evidence_df = pd.read_csv(mq_folder_path / "evidence.txt", sep="\t")

    # Remove contaminants and decoys
    evidence_df = evidence_df[
        (evidence_df["Reverse"] != "+") & (evidence_df["Potential contaminant"] != "+")
    ]

    evidence_df["Modified sequence"] = evidence_df["Modified sequence"].str.replace(
        "(Acetyl (Protein N-term))", "M", regex=False
    )
    evidence_df["Modified sequence"] = evidence_df["Modified sequence"].str.replace(
        "(Oxidation (M))", "M", regex=False
    )

    evidence_df["Modified sequence"] = evidence_df["Modified sequence"].str.replace(
        r"(S|T|Y)(Phospho (STY))", lambda x: "p" + x.group(1), regex=True
    )
    return evidence_df


def plot_evidence_qc_plots(
    evidence_df: pd.DataFrame, output_path: Path, batch_name: str, prefix: str
):
    plot_miscleavages(evidence_df, output_path, batch_name)

    peptide_intensities_per_channel = get_intensity_per_channel(evidence_df)
    plot_evidence_qc_plots_grouped(
        peptide_intensities_per_channel,
        output_path,
        batch_name,
        prefix,
        xlabel="Channel",
    )

    peptide_intensities_per_fraction = get_intensity_per_fraction(evidence_df)
    plot_evidence_qc_plots_grouped(
        peptide_intensities_per_fraction,
        output_path,
        batch_name,
        prefix,
        xlabel="Fraction",
    )

    return peptide_intensities_per_fraction, peptide_intensities_per_channel


def get_intensity_per_channel(evidence_df: pd.DataFrame):
    # sum reporter intensity of same modified sequence per batch
    df = evidence_df.set_index("Modified sequence")
    df = df.filter(like="Reporter intensity corrected")
    df = df.groupby("Modified sequence").agg("sum")

    df = df.fillna(0)
    df = df.loc[df.sum(axis=1) > 0]
    df = df.replace(0.0, np.nan)

    df = df.rename(columns=lambda x: x.replace("Reporter intensity corrected ", ""))

    return df


def get_intensity_per_fraction(evidence_df: pd.DataFrame):
    # sum reporter intensity of same modified sequence per batch
    df = evidence_df[["Modified sequence", "Fraction", "Intensity"]]
    df = df.pivot_table(
        index="Modified sequence", columns="Fraction", values="Intensity", aggfunc="sum"
    )

    df = df.fillna(0)
    df = df.loc[df.sum(axis=1) > 0]
    df = df.replace(0.0, np.nan)

    return df


def plot_summary_by_fraction(
    summary_df: pd.DataFrame,
    output_path: Path,
    column: str,
    ylim: float,
    title_suffix: str = "",
):
    plt.figure(figsize=(8, 6))
    ax = sns.barplot(
        data=summary_df,
        x="Fraction",
        y=column,
        hue="Fraction",
        palette="colorblind",
        legend=False,
    )
    plt.ylim(0, ylim)
    ax.tick_params(axis="x", rotation=90)
    plt.title(
        f"{summary_df.iloc[0]['Experiment']} - Average {round(summary_df[column].mean(), 0)}{title_suffix}"
    )
    plt.tight_layout()
    plt.savefig(output_path / f"{column.replace(' ', '_').replace('/', '')}.png")
    plt.close()


def plot_miscleavages(evidence_df: pd.DataFrame, output_path: Path, batch_name: str):
    miscleavages = evidence_df["Missed cleavages"].value_counts()

    plt.figure()
    ax = sns.barplot(x=miscleavages.index, y=miscleavages.values)
    ax.set_xticks(range(miscleavages.index.max() + 1), labels=miscleavages.index)
    plt.ylabel("Counts")
    plt.xlabel("Miscleavages")
    plt.title(batch_name)
    plt.savefig(output_path / "Miscleavages.png")
    plt.close()


def plot_evidence_qc_plots_grouped(
    intensities_per_channel: pd.DataFrame,
    output_path: Path,
    batch_name: str,
    prefix: str,
    xlabel: str,
):
    num_ids_per_channel = intensities_per_channel.count()
    summed_intensity_per_channel = intensities_per_channel.sum()

    plt.figure(figsize=(8, 6))
    plt.bar(num_ids_per_channel.index, num_ids_per_channel.values)
    plt.ylabel("IDs")
    plt.xlabel(xlabel)
    plt.title(batch_name)
    plt.tight_layout()
    plt.savefig(output_path / f"{prefix}_ID_{xlabel}.png")
    plt.close()

    plt.figure(figsize=(8, 6))
    plt.bar(summed_intensity_per_channel.index, summed_intensity_per_channel.values)
    plt.ylabel("Summed intensity")
    plt.xlabel(xlabel)
    plt.title(batch_name)
    plt.tight_layout()
    plt.savefig(output_path / f"{prefix}_sumInt_{xlabel}.png")
    plt.close()

    plt.figure(figsize=(8, 6))
    ax = sns.boxplot(data=np.log2(intensities_per_channel))
    ax.tick_params(axis="x", rotation=90)
    plt.ylabel("log2(Intensity)")
    plt.xlabel(xlabel)
    plt.title(batch_name)
    plt.tight_layout()
    plt.savefig(output_path / f"{prefix}_boxInt_{xlabel}.png")
    plt.close()


def plot_imac_selectivity_per_fraction(
    evidence_df: pd.DataFrame,
    evidence_df_only_phosho: pd.DataFrame,
    output_path: Path,
    batch_name: str,
):
    imac_selectivity_df = get_imac_selectivity(evidence_df, evidence_df_only_phosho)

    plt.figure(figsize=(8, 6))
    ax = sns.barplot(
        x="Fraction",
        hue="Fraction",
        y="Selectivity IMAC",
        data=imac_selectivity_df,
        palette="colorblind",
        legend=False,
    )
    ax.tick_params(axis="x", rotation=90)
    plt.title(batch_name)
    plt.ylim(0, 1)
    plt.axhline(y=0.8, color="red", linestyle="--")
    plt.tight_layout()
    plt.savefig(output_path / "IMAC_selectivity.png")
    plt.close()

    return imac_selectivity_df


def get_imac_selectivity(
    evidence_df: pd.DataFrame,
    evidence_df_only_phosho: pd.DataFrame,
):
    intensity_per_fraction_all = get_intensity_per_fraction(evidence_df)
    intensity_per_fraction_only_phospho = get_intensity_per_fraction(
        evidence_df_only_phosho
    )

    imac_selectivity_df = (
        intensity_per_fraction_only_phospho.sum() / intensity_per_fraction_all.sum()
    )
    imac_selectivity_df = imac_selectivity_df.to_frame().reset_index()
    imac_selectivity_df.columns = ["Fraction", "Selectivity IMAC"]

    return imac_selectivity_df


def read_protein_groups_txt(mq_folder_path: Path):
    protein_groups_df = pd.read_csv(mq_folder_path / "proteinGroups.txt", sep="\t")

    # Remove contaminants and decoys
    protein_groups_df = protein_groups_df[
        (protein_groups_df["Reverse"] != "+")
        & (protein_groups_df["Potential contaminant"] != "+")
    ]

    protein_groups_df = protein_groups_df[protein_groups_df["Gene names"].notnull()]
    protein_groups_df = protein_groups_df.set_index("Gene names")

    protein_groups_df = protein_groups_df.filter(
        like="Reporter intensity corrected"
    ).dropna(how="all")
    protein_groups_df = protein_groups_df.rename(
        columns=lambda x: x.replace("Reporter intensity corrected ", "").split()[0]
    )

    protein_groups_df = protein_groups_df.replace(0, np.nan)

    return protein_groups_df


def group_evidence_by_gene_name(evidence_df: pd.DataFrame):
    protein_groups_df = evidence_df[evidence_df["Gene names"].notnull()]
    protein_groups_df = protein_groups_df.groupby("Gene names").agg(
        "sum", numeric_only=True
    )

    protein_groups_df = protein_groups_df.replace(0, np.nan)

    protein_groups_df = protein_groups_df.filter(
        like="Reporter intensity corrected"
    ).dropna(how="all")
    protein_groups_df = protein_groups_df.rename(
        columns=lambda x: x.replace("Reporter intensity corrected ", "").split()[0]
    )

    return protein_groups_df


def read_basket_annotation_file(basket_table_path: Path):
    return pd.read_csv(
        basket_table_path,
        sep="\t",
        na_values="",
    )


def plot_basket_identification_rates(
    basket_members: pd.DataFrame, protein_groups_df: pd.DataFrame, output_path: Path
):
    basket_sizes = basket_members.notnull().sum()

    plt.figure(figsize=(10, 10))
    for i, (col, basket_size) in enumerate(zip(basket_members.columns, basket_sizes)):
        plt.subplot(3, 3, i + 1)

        basket_members_detected = protein_groups_df[
            protein_groups_df.index.isin(basket_members[col].dropna())
        ]
        detection_rate = basket_members_detected.count().astype(float) / basket_size

        ax = sns.barplot(
            x=detection_rate.index,
            hue=detection_rate.index,
            y=detection_rate.values,
            palette="colorblind",
            legend=False,
        )

        # ax.bar_label(ax.containers[0], labels=basket_members_detected.count())

        basket_counts = (
            basket_members_detected.count()
        )  # Ensure this is a list-like object

        # Loop over each bar container and apply the labels
        for i, container in enumerate(ax.containers):
            # Ensure there are enough labels for each container
            if len(basket_counts) > i:
                ax.bar_label(
                    container, labels=[basket_counts.iloc[i]]
                )  # Apply the correct label for each container
            else:
                logger.info(f"Not enough labels for container {i}")

        plt.ylim(0, 1)
        plt.ylabel("% basket identified")
        plt.xlabel("Channel")
        plt.title(col)
    plt.tight_layout()
    plt.savefig(output_path / "basket_IDs.png")
    plt.close()


def write_batch_statistics_table(
    batch_name: str,
    peptide_intensities_per_fraction: pd.DataFrame,
    peptide_intensities_per_channel: pd.DataFrame,
    ppeptide_intensities_per_fraction: pd.DataFrame,
    ppeptide_intensities_per_channel: pd.DataFrame,
    imac_selectivity_df: pd.DataFrame,
    protein_groups_df: pd.DataFrame,
    summary_df_totals: pd.Series,
    summary_df: pd.DataFrame,
    output_path: Path,
    search_type: str,
):
    # stats per fraction
    stats_by_fraction_summary_df = summary_df.set_index("Fraction")[
        ["PSMs", "ID rate", "MS scans", "MS/MS scans"]
    ]

    peptide_id_count = peptide_intensities_per_fraction.count().tolist()
    peptide_intensity_sum = peptide_intensities_per_fraction.sum().tolist()

    ppeptide_stats = {}
    if search_type == "PP":
        ppeptide_id_count = ppeptide_intensities_per_fraction.count().tolist()
        ppeptide_intensity_sum = ppeptide_intensities_per_fraction.sum().tolist()
        imac_selectivity = (imac_selectivity_df["Selectivity IMAC"] * 100).tolist()

        ppeptide_stats = {
            "Mod_phosphopeptides": ppeptide_id_count,
            "Summed phosphopeptide intensity": ppeptide_intensity_sum,
            "IMAC selectivity": imac_selectivity,
        }

    stats_by_fraction_df = pd.DataFrame(
        {
            "Mod_peptides": peptide_id_count,
            "Summed peptide intensity": peptide_intensity_sum,
        }
        | ppeptide_stats,
        index=peptide_intensities_per_fraction.columns,
    )
    # make indices strings
    stats_by_fraction_df.index = stats_by_fraction_df.index.astype(str)
    stats_by_fraction_summary_df.index = stats_by_fraction_summary_df.index.astype(str)

    # we had a case of a fraction split into 2 raw files as machine stopped:
    if stats_by_fraction_summary_df.index.tolist() != stats_by_fraction_df.index.tolist():
        stats_by_fraction_summary_df = add_postfix_to_fraction_duplicates(stats_by_fraction_summary_df)

    stats_by_fraction_df = pd.merge(
        stats_by_fraction_summary_df,
        stats_by_fraction_df,
        left_index=True,
        right_index=True,
        how="outer"
    )
    stats_by_fraction_df = stats_by_fraction_df.reset_index()
    stats_by_fraction_df.insert(1, "Channel", np.nan)

    # stats per channel
    protein_intensity_sum = protein_groups_df.notnull().sum()
    peptide_id_count = peptide_intensities_per_channel.count().tolist()
    peptide_intensity_sum = peptide_intensities_per_channel.sum().tolist()

    ppeptide_stats = {}
    if search_type == "PP":
        ppeptide_id_count = ppeptide_intensities_per_channel.count().tolist()
        ppeptide_intensity_sum = ppeptide_intensities_per_channel.sum().tolist()

        ppeptide_stats = {
            "Mod_phosphopeptides": ppeptide_id_count,
            "Summed phosphopeptide intensity": ppeptide_intensity_sum,
        }

    stats_by_channel_df = pd.DataFrame(
        {
            "Channel": peptide_intensities_per_channel.columns.tolist(),
            "Proteins": protein_intensity_sum,
            "Mod_peptides": peptide_id_count,
            "Summed peptide intensity": peptide_intensity_sum,
        }
        | ppeptide_stats,
        index=peptide_intensities_per_channel.columns.tolist(),
    )

    # stats total
    ppeptide_stats = {}
    if search_type == "PP":
        ppeptide_stats = {
            "Mod_phosphopeptides": [len(ppeptide_intensities_per_fraction)],
            "Summed phosphopeptide intensity": [
                ppeptide_intensities_per_channel.sum().sum()
            ],
            "IMAC selectivity": [imac_selectivity_df["Selectivity IMAC"].mean() * 100],
        }

    totals_df = pd.DataFrame(
        {
            "Channel": ["total"],
            "Fraction": ["total"],
            "Proteins": [len(protein_groups_df)],
            "Mod_peptides": [len(peptide_intensities_per_fraction)],
            "PSMs": [summary_df_totals["PSMs"]],
            "ID rate": [summary_df["ID rate"].mean()],
            "MS scans": [summary_df_totals["MS scans"]],
            "MS/MS scans": [summary_df_totals["MS/MS scans"]],
            "Summed peptide intensity": [peptide_intensities_per_channel.sum().sum()],
        }
        | ppeptide_stats,
        index=["total"],
    )

    table = pd.concat([stats_by_fraction_df, stats_by_channel_df, totals_df])
    table["Experiment"] = batch_name

    table.to_csv(
        output_path / f"ID_table_{search_type}.csv",
        index=False,
        float_format="%.1f",
    )


def add_postfix_to_fraction_duplicates(stats_by_fraction_summary_df: pd.DataFrame):
    # find the duplicated indices
    duplicated_indices = stats_by_fraction_summary_df.index[
        stats_by_fraction_summary_df.index.duplicated()
    ].unique()
    logger.warning(
        f"Duplication in fraction indices found: {duplicated_indices.tolist()}."
    )
    # add postfix to duplicated indices
    for dup_index in duplicated_indices:
        dup_count = 0
        for i in range(len(stats_by_fraction_summary_df.index)):
            if stats_by_fraction_summary_df.index[i] == dup_index:
                dup_count += 1
                stats_by_fraction_summary_df.index.values[i] = (
                    f"{dup_index}_{dup_count}"
                )
    return stats_by_fraction_summary_df


def append_channel_statistics_to_tracking_table(
    batch_name: str,
    peptide_intensities_per_channel: pd.DataFrame,
    protein_numbers: pd.DataFrame,
    output_path: Path,
    search_type: str,
):

    prefix = "pept"
    if search_type == "PP":
        prefix = "ppept"

    # Create DataFrame for IDs
    IDs = pd.DataFrame(
        {
            "Batch": [batch_name] * 11,
            "Channel": np.arange(1, 12),
            f"{prefix}_id": peptide_intensities_per_channel.count(),
            f"{prefix}_int_mean": peptide_intensities_per_channel.mean().round(1),
            f"{prefix}_int": peptide_intensities_per_channel.sum().round(1),
            f"{prefix}_proteins": protein_numbers.count(),
        }
    )
    tracking_table_path = output_path / f"combined_{search_type}_ID.csv"
    if not tracking_table_path.is_file():
        IDs.to_csv(
            tracking_table_path,
            index=False,
            float_format="%.1f",
        )
        return IDs
    else:
        peptide_intensities_per_channel = pd.read_csv(
            tracking_table_path,
            sep=",",
        )
        peptide_intensities_per_channel = pd.concat(
            [peptide_intensities_per_channel, IDs], ignore_index=True
        )
        peptide_intensities_per_channel = peptide_intensities_per_channel.drop_duplicates()
        peptide_intensities_per_channel.to_csv(
            tracking_table_path,
            index=False,
            float_format="%.1f",
        )
    return peptide_intensities_per_channel


def plot_tracking_table_pdf(
    tracking_table: pd.DataFrame, output_path: Path, search_type: str
):
    prefix = "pept"
    if search_type == "PP":
        prefix = "ppept"

    _, axes = plt.subplots(2, 1, figsize=(15, 10))

    sns.scatterplot(
        ax=axes[0],
        data=tracking_table,
        x=tracking_table.index,
        y=f"{prefix}_id",
        hue="Batch",
        legend=False,
        palette="tab10",
        s=10,
    )
    sns.scatterplot(
        ax=axes[1],
        data=tracking_table,
        x=tracking_table.index,
        y=f"{prefix}_int_mean",
        hue="Batch",
        legend=False,
        palette="tab10",
        s=10,
    )
    plt.yscale("log")

    xticks = np.rint(np.linspace(0, len(tracking_table.index) - 1, 50))
    labels = tracking_table["Batch"] + "_" + tracking_table["Channel"].astype(str)
    axes[1].set_xticks(xticks, labels=labels.iloc[xticks], rotation=90)

    plt.tight_layout()
    plt.savefig(output_path / f"combined_{search_type}.pdf")
    plt.close()


if __name__ == "__main__":
    mq_folder_base_path = sys.argv[1]
    search_type = sys.argv[2]  # FP or PP
    main(mq_folder_base_path, search_type)
