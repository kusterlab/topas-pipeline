from pathlib import Path

import pandas as pd
 
from topas_pipeline import picked_group, config, sample_annotation
from topas_pipeline import preprocess as prep

CONFIG_FILE_PATH = './data/test_config.json'


def test_run_picked_group_fdr():
    """Check if results from picked_group_fdr.main() are the same as in the reference run.

    Takes approximately 50 seconds.
    """
    configs = config.load(CONFIG_FILE_PATH)

    results_folder = configs["results_folder"]
    fasta_file = configs["preprocessing"]["fasta_file"]
    evidence_file = Path(results_folder) / "evidence.txt"
    picked_output_file = (
        Path(results_folder) / "integration_tests" / "pickedGeneGroups.txt"
    )
    if picked_output_file.is_file():
        picked_output_file.unlink()

    picked_group.run_picked_group_fdr(
        str(evidence_file), str(picked_output_file), fasta_file
    )

    results_df = pd.read_csv(picked_output_file, sep="\t")
    results_df = results_df.sort_values(by="Protein IDs")
    results_df = results_df[results_df["Q-value"] < 0.01]

    output_file_ref = (
        Path(results_folder) / "integration_tests" / "pickedGeneGroups_ref.txt"
    )
    results_ref_df = pd.read_csv(output_file_ref, sep="\t")
    results_ref_df = results_ref_df.sort_values(by="Protein IDs")

    pd.testing.assert_frame_equal(results_df, results_ref_df, atol=1e-5)


def test_do_quant():
    """Check if results from picked_group_fdr.do_quant are the same as in the reference run.

    Takes approximately 150 seconds.
    """
    configs = config.load(CONFIG_FILE_PATH)
    results_folder = configs["results_folder"]
    fasta_file = configs["preprocessing"]["fasta_file"]    

    picked_output_file_with_quant = (
        Path(results_folder) / "integration_tests" / "pickedGeneGroups_with_quant.txt"
    )
    if picked_output_file_with_quant.is_file():
        picked_output_file_with_quant.unlink()

    sample_annotation_df = sample_annotation.load_sample_annotation(
        configs["sample_annotation"]
    )

    df = prep.load_sample_data(
        configs["results_folder"],
        sample_annotation_df,
        configs["simsi"]["simsi_folder"],
        configs["preprocessing"]["raw_data_location"],
        configs["preprocessing"]["run_simsi"],
        configs["preprocessing"]["run_lfq"],
        configs["preprocessing"]["debug"],
        "fp",
    )

    picked_output_file_ref = (
        Path(results_folder) / "integration_tests" / "pickedGeneGroups_ref.txt"
    )
    num_threads = 4
    picked_group.do_quant(
        df,
        str(picked_output_file_ref),
        fasta_file,
        str(picked_output_file_with_quant),
        num_threads,
    )

    results_with_quant_df = pd.read_csv(picked_output_file_with_quant, sep="\t")

    picked_output_file_with_quant_ref = (
        Path(results_folder)
        / "integration_tests"
        / "pickedGeneGroups_with_quant_ref.txt"
    )
    results_with_quant_ref_df = pd.read_csv(picked_output_file_with_quant_ref, sep="\t")

    pd.testing.assert_frame_equal(results_with_quant_df, results_with_quant_ref_df)


def test_remap_gene_names():
    configs = config.load(CONFIG_FILE_PATH)
    fasta_file = configs["preprocessing"]["fasta_file"]

    sample_annotation_df = sample_annotation.load_sample_annotation(
        configs["sample_annotation"]
    )

    df = prep.load_sample_data(
        configs["results_folder"],
        sample_annotation_df,
        configs["simsi"]["simsi_folder"],
        configs["preprocessing"]["raw_data_location"],
        configs["preprocessing"]["run_simsi"],
        configs["preprocessing"]["run_lfq"],
        configs["preprocessing"]["debug"],
        "pp",
    )

    genes_before = df["Gene names"]

    assert "WARS" in genes_before.values
    assert "WARS1" not in genes_before.values

    remapped_genes_df = picked_group.remap_gene_names(df, fasta_file)

    genes_after = remapped_genes_df["Gene names"]

    assert len(genes_after) == len(genes_before)

    assert "WARS" not in genes_after.values
    assert "WARS1" in genes_after.values


if __name__ == "__main__":
    test_do_quant()
