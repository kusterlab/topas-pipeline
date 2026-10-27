import argparse
import logging
import pandas as pd
import os
import json

logger = logging.getLogger("ScoreCompare")
logging.basicConfig(level=logging.INFO)


def get_poi_df(path_annot_fp, path_z_measures_fp, poi_col='POI_REPORT', id_col='Gene names'):
    df_annot_fp = pd.read_csv(path_annot_fp)
    df_z_measures_fp = pd.read_csv(path_z_measures_fp, sep="\t")
    poi_genes_fp = df_annot_fp[df_annot_fp[poi_col].notna()][id_col].unique().tolist()
    return df_z_measures_fp[df_z_measures_fp[id_col].isin(poi_genes_fp)].set_index(id_col).T.reset_index().rename(columns={'index': id_col})

def generate_plot_data_kinase_score_change(df_before, df_after, id_col="Sample name", threshold=1.0):
    """
    Generate a DataFrame containing the change in kinase scores between two conditions.

    Parameters:
    - df_before: DataFrame containing kinase scores before the condition.
    - df_after: DataFrame containing kinase scores after the condition.
    - id_col: Column name that uniquely identifies each sample (default is "Sample name").
    - threshold: Minimum absolute score to consider a kinase score as significant (default is 1.0).

    Returns:
    - plot_data: DataFrame containing the change in kinase scores for significant kinases.
    - rank_order: List of kinases ordered by their maximum absolute change.
    """
    
    # Treat every column except the sample ID as a kinase score.
    kinase_cols = [col for col in df_before.columns if col != id_col]
    after_kinase_cols = [col for col in df_after.columns if col != id_col]

    if set(kinase_cols) != set(after_kinase_cols):
        raise ValueError("Before and after files must contain the same kinase columns.")

    # Ensure scores are numeric.
    df_before[kinase_cols] = df_before[kinase_cols].apply(pd.to_numeric, errors="coerce")
    df_after[kinase_cols] = df_after[kinase_cols].apply(pd.to_numeric, errors="coerce")

    # Pair measurements by sample. This also checks for duplicate or unmatched IDs.
    paired = df_before[[id_col] + kinase_cols].merge(
        df_after[[id_col] + kinase_cols],
        on=id_col,
        how="inner",
        suffixes=("_before", "_after"),
    )
    paired = paired.fillna(0.00)
    paired.set_index(id_col, inplace=True)

    # Positive ΔZ means the score increased after; negative means it decreased.
    delta = pd.DataFrame(
        {
            kinase: paired[f"{kinase}_after"] - paired[f"{kinase}_before"]
            for kinase in kinase_cols
        },
        index=paired.index,
    )

    # Mark a sample/kinase pair when either score reaches the cutoff.
    important = pd.DataFrame(
        {
            kinase: (
                paired[f"{kinase}_before"].abs().ge(threshold)
                | paired[f"{kinase}_after"].abs().ge(threshold)
            )
            for kinase in kinase_cols
        },
        index=paired.index,
    )

    # Keep kinases that meet the cutoff in at least one paired sample.
    important_kinases = important.columns[important.any(axis=0)].tolist()

    if not important_kinases:
        raise ValueError(f"No kinases met the |Z| >= {threshold} cutoff.")

    # Mask patient/kinase pairs that did not meet the cutoff.
    plot_data = delta[important_kinases].where(important[important_kinases])

    # Order kinase columns by median absolute change.
    rank_order = (
        plot_data.abs()
        .max()
        .sort_values(ascending=False)
        .index
        .tolist()
    )
    summary_df = (
        plot_data[rank_order]
        .rename_axis(index=id_col, columns="Kinase")
        .reset_index()
        .melt(id_vars=id_col, var_name="Kinase", value_name="Delta Z")
        .dropna()
    )

    return summary_df, rank_order


def generate_strip_plot(summary, rank_order, title, id_col, xlabel, ylabel="Change in Z score (after - before)"):
    """
    Generate a strip plot of kinase score changes.

    Parameters:
    - summary: DataFrame containing the change in kinase scores.
    - rank_order: List of kinases ordered by their maximum absolute change.
    - id_col: Column name that uniquely identifies each sample.

    Returns:
    - strip_fig: Plotly figure object representing the strip plot.
    """
    import plotly.express as px

    # Strip plot: each point is one patient's ΔZ
    strip_fig = px.strip(
        summary,
        x="Kinase",
        y="Delta Z",
        category_orders={"Kinase": rank_order},
        hover_data=[id_col],
        title=title,
        labels={"Delta Z": ylabel, "Kinase": xlabel},
    )
    strip_fig.update_traces(jitter=0.25, marker={"size": 7, "opacity": 0.75})
    strip_fig.update_layout(xaxis_tickangle=-60)

    return strip_fig

def generate_box_plot(summary, rank_order, title, id_col, xlabel, ylabel="Change in Z score (after - before)"):
    import plotly.express as px

    box_fig = px.box(
        summary,
        x="Kinase",
        y="Delta Z",
        category_orders={"Kinase": rank_order},
        hover_data=[id_col],
        title=title,
        labels={"Delta Z": ylabel, "Kinase": xlabel},
    )
    box_fig.update_layout(xaxis_tickangle=-60)

    return box_fig

def plot_change_ck_scores(baseline_results_dir, current_results_dir, threshold=1.0, plot_type="strip"):
    ck_scores_file_path_before = os.path.join(baseline_results_dir, "topas_scores", "ck_substrate_phosphorylation_scores_expressioncorrected.tsv")
    ck_scores_file_path_after = os.path.join(current_results_dir, "topas_scores", "ck_substrate_phosphorylation_scores_expressioncorrected.tsv")
    df_before = pd.read_csv(ck_scores_file_path_before, sep="\t")
    df_after = pd.read_csv(ck_scores_file_path_after, sep="\t")
    summary_df, rank_order = generate_plot_data_kinase_score_change(df_before, df_after, id_col="Sample name", threshold=threshold)
    fig = PLOT_TYPE_REGISTRY[plot_type](summary_df, rank_order, title="Patient-level CK score changes", id_col="Sample name", xlabel="Kinase")
    return fig

def plot_change_rtk_scores(baseline_results_dir, current_results_dir, threshold=1.0, plot_type="strip"):
    rtk_scores_file_path_before = os.path.join(baseline_results_dir, "topas_scores", "topas_rtk_scores_zscored.tsv")
    rtk_scores_file_path_after = os.path.join(current_results_dir, "topas_scores", "topas_rtk_scores_zscored.tsv")
    df_before = pd.read_csv(rtk_scores_file_path_before, sep="\t")
    df_after = pd.read_csv(rtk_scores_file_path_after, sep="\t")
    summary_df, rank_order = generate_plot_data_kinase_score_change(df_before, df_after, id_col="Sample", threshold=threshold)
    fig = PLOT_TYPE_REGISTRY[plot_type](summary_df, rank_order, title="Patient-level RTK score changes", id_col="Sample", xlabel="Kinase")
    return fig

def plot_change_poi_scores(baseline_results_dir, current_results_dir, threshold=1.0, plot_type="strip"):
    path_annot_fp_before = os.path.join(baseline_results_dir, "annot_fp.csv")
    path_z_measures_fp_before = os.path.join(baseline_results_dir, "full_proteome_measures_z.tsv")
    path_annot_fp_after = os.path.join(current_results_dir, "annot_fp.csv")
    path_z_measures_fp_after = os.path.join(current_results_dir, "full_proteome_measures_z.tsv")
    poi_col = 'POI_REPORT'
    df_poi_before = get_poi_df(path_annot_fp_before, path_z_measures_fp_before, poi_col=poi_col)
    df_poi_after = get_poi_df(path_annot_fp_after, path_z_measures_fp_after, poi_col=poi_col)
    summary_df, rank_order = generate_plot_data_kinase_score_change(df_poi_before, df_poi_after, id_col="Gene names", threshold=threshold)
    fig = PLOT_TYPE_REGISTRY[plot_type](summary_df, rank_order, title=f"Patient-level {poi_col} score changes", id_col="Gene names", xlabel="Gene names")
    return fig

PLOT_TYPE_REGISTRY = {
    "strip": generate_strip_plot,
    "box": generate_box_plot,
}

COMPARISON_PLOT_REGISTRY = {
    "ck": plot_change_ck_scores,
    "rtk": plot_change_rtk_scores,
    "poi": plot_change_poi_scores,
}

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--baseline", required=True, type=str, help="Path to the baseline Topas pipeline results directory"
    )
    parser.add_argument(
        "--current", required=True, type=str, help="Path to the current Topas pipeline results directory"
    )
    parser.add_argument(
        "--output", required=True, type=str, help="Path to the output directory to save the generated strip plot HTML file"
    )
    parser.add_argument(
        "--threshold", type=float, default=1.0, help="Threshold for significant kinase score changes."
    )
    args = parser.parse_args()
    return args

def main(args):
    logger.info(f"Configuration for score comparison:\n{json.dumps(COMPARISON_CONFIG, indent=4)}")
    logger.info("Starting score comparison plots...")
    plot_config = COMPARISON_CONFIG["plots"]
    os.makedirs(args.output, exist_ok=True)
    for key, plot_type in plot_config.items():
        logger.info(f"Generating {key} score change plot with {plot_type} plot type...")
        fig = COMPARISON_PLOT_REGISTRY[key](args.baseline, args.current, threshold=args.threshold, plot_type=plot_type)
        fig.write_html(os.path.join(args.output, f"{key}_score_changes.html"))
        logger.info(f"Saved {key} score change plot to {os.path.join(args.output, f'{key}_score_changes.html')}")
    logger.info("Completed all score comparison plots.")

if __name__ == "__main__":
    
    args = parse_args()
    COMPARISON_CONFIG = {
        "baseline_results_dir": args.baseline,
        "current_results_dir": args.current,
        "threshold": args.threshold,
        "output": args.output,
        "plots": {
            "ck": "box",
            "rtk": "box",
            "poi": "box"
        }
    }
    main(args)
    