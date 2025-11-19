import sys
import pandas as pd

# python tests/integration_tests/compare_tables.py "results/2025.04.01_example_run" "results/2025.04.01_example_run2"

PIPELINE_OUTPUT_FILES = [
    "preprocessed_fp.csv",
    # "preprocessed_pp.csv",
    "annot_fp.csv",
    # "annot_pp.csv",
    "full_proteome_measures_rank.tsv",
    "full_proteome_measures_fc.tsv",
    "full_proteome_measures_z.tsv",
    # "phospho_measures_rank.tsv",
    # "phospho_measures_fc.tsv",
    # "phospho_measures_z.tsv",
    "topas_scores/ck_substrate_phosphorylation_scores_expressioncorrected.tsv",
    "topas_scores/rtk_substrate_phosphorylation_scores.tsv",
    "topas_scores/protein_phosphorylation_scores.tsv",
    "topas_scores/topas_rtk_scores_zscored.tsv",
]


def main(argv):
    folder1 = argv[0]
    folder2 = argv[1]
    for file_name in PIPELINE_OUTPUT_FILES:
        print(f"Comparing {file_name}")
        sep = "\t"
        if file_name.endswith(".csv"):
            sep = ","

        df1 = pd.read_csv(f"{folder1}/{file_name}", sep=sep)
        df2 = pd.read_csv(f"{folder2}/{file_name}", sep=sep)

        pd.testing.assert_frame_equal(df1, df2, rtol=1e-3, atol=1e-2, check_dtype=False)


if __name__ == "__main__":
    main(sys.argv[1:])