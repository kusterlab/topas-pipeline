import pandas as pd

# python tests/integration_tests/compare_tables.py

folder1 = "results/2025.04.01_example_run"
folder2 = "results/2025.04.01_example_run_2"

file_names = [
    "basket_scores_4th_gen.tsv",
    "kinase_results/kinase_scores.tsv",
    "protein_results/protein_scores.tsv",
    "preprocessed_fp.csv",
    "preprocessed_pp.csv",
    "annot_fp.csv",
    "annot_pp.csv",
    "full_proteome_measures_rank.tsv",
    "full_proteome_measures_fc.tsv",
    "full_proteome_measures_z.tsv",
    "full_proteome_measures_p.tsv",
    "phospho_measures_rank.tsv",
    "phospho_measures_fc.tsv",
    "phospho_measures_z.tsv",
    "phospho_measures_p.tsv",
]

for file_name in file_names:
    print(f"Comparing {file_name}")
    sep = "\t"
    if file_name.endswith(".csv"):
        sep = ","

    df1 = pd.read_csv(f"{folder1}/{file_name}", sep=sep)
    df2 = pd.read_csv(f"{folder2}/{file_name}", sep=sep)

    pd.testing.assert_frame_equal(df1, df2, rtol=1e-3, atol=1e-2, check_dtype=False)
