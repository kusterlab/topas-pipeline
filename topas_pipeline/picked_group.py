import sys
import re
import os
import json
import logging
from pathlib import Path
from typing import List

import pandas as pd
import numpy as np
import argparse

from . import config
from . import utils
from . import sample_annotation
import topas_pipeline.preprocess as pre

import picked_group_fdr.picked_group_fdr as picked
import picked_group_fdr.quantification as quant
import picked_group_fdr.digest as digest
import picked_group_fdr.protein_annotation as protein_annotation
import picked_group_fdr.helpers as helpers
import picked_group_fdr.parsers.maxquant as maxquant
import picked_group_fdr.columns as columns
from picked_group_fdr.results import ProteinGroupResults
from picked_group_fdr.protein_groups import ProteinGroups
from picked_group_fdr.precursor_quant import PrecursorQuant

# hacky way to get the package logger instead of just __main__ when running as a module
logger = logging.getLogger(__package__ + "." + Path(__file__).stem)


def prepare_evidence_for_picked(df, evidence_file):
    """
    Columns needed for picked protein group FDR:
    - Modified sequence
    - Leading proteins (not used for all methods, but should still be in the input file)
    - PEP
    - Experiment
    - Charge
    - Intensity
    - Raw file
    - Fraction
    - Id
    """
    logger.info("Generating evidence file for picked group FDR")
    info_cols = [
        "Modified sequence",
        "Proteins",
        "Leading proteins",
        "Gene names",
        "Raw file",
        "Fraction",
        "Batch",
        "Charge",
        "PEP",
        "id",
    ]
    if "Gene names" not in df.columns:
        info_cols.remove("Gene names")
    intensity_cols = ["Intensity"] + df.columns[
        df.columns.str.startswith("Reporter intensity")
    ].tolist()

    df = df[info_cols + intensity_cols]

    df.fillna({"PEP": "NaN"}).to_csv(evidence_file, sep="\t", index=False)


def run_picked_group_fdr(
    evidence_file, picked_fdr_file, fasta_file, gene_level: bool = True
):
    if os.path.exists(picked_fdr_file):
        logger.info(
            f"Reusing previously generated results from picked protein group FDR"
        )
        return

    logger.info("Running picked group FDR")

    picked.main(
        [
            "--mq_evidence",
            evidence_file,
            "--protein_groups_out",
            picked_fdr_file,
            "--fasta",
            fasta_file,
            "--enzyme",
            "trypsinp",
            "--min-length",
            "6",
            "--cleavages",
            "3",
            "--lfq_min_peptide_ratios",
            "1",
            "--gene_level",
            "--suppress_missing_peptide_warning",
            "--methods",
            "picked_protein_group_mq_input",
        ]
    )


def picked_protein_grouping(
    df: pd.DataFrame,
    results_folder: str,
    protein_fdr_cutoff: float,
    fasta_file: str,
    fdr_num_threads: int,
    evidence_file_base: str = "evidence",
    output_file_base: str = "pickedGeneGroups",
):
    # propagate logs from picked_group_fdr package to our current log handlers
    logging.getLogger("picked_group_fdr").handlers = logging.getLogger(
        __package__
    ).handlers

    # Picked protein grouping
    evidence_file = f"{results_folder}/{evidence_file_base}.txt"

    if os.path.exists(evidence_file):
        logger.info(
            f"Reusing previously generated evidence file for picked protein group FDR"
        )
    else:
        prepare_evidence_for_picked(df, evidence_file)

    picked_fdr_file = f"{results_folder}/{output_file_base}.txt"
    gene_level = True
    if "Gene names" not in df.columns:
        gene_level = False
    run_picked_group_fdr(evidence_file, picked_fdr_file, fasta_file, gene_level)
    picked_fdr_file_with_quant = f"{results_folder}/{output_file_base}_with_quant.txt"
    do_quant(
        df, picked_fdr_file, fasta_file, picked_fdr_file_with_quant, fdr_num_threads
    )

    # read in quantified genes
    df = pd.read_csv(picked_fdr_file_with_quant, sep="\t")

    # Filter at 1% gene-level FDR
    before = df.shape[0]
    df = df[df["Q-value"] < protein_fdr_cutoff]
    after = df.shape[0]

    log_file = open(results_folder + "/log_file.txt", "a")
    log_file.write("\n Full proteome FDR filtering (Q-value): Before, after\n")
    log_file.write(str(before) + ", " + str(after) + "\n")
    log_file.close()

    df = df.rename(
        columns=lambda x: x.replace(
            "LFQ Intensity Reporter intensity", "Reporter intensity"
        )
    )

    # only report genes if they have >=50% of the identified peptides for the gene group
    # (a.k.a. Majority protein). An example where this makes sense is AXIN1, which
    # is sometimes grouped with MACF1 but only has ~2 peptides, whereas MACF1 has ~3000.
    df["Gene names"] = df["Majority protein IDs"]

    return df


def do_quant(
    df, picked_fdr_file, fasta_file, picked_fdr_file_with_quant, fdr_num_threads: int
):
    if os.path.exists(picked_fdr_file_with_quant):
        logger.info(
            f"Reusing previously generated results from picked protein group FDR with quantification"
        )
        return

    args = quant.parse_args(
        [
            "--mq_evidence",
            "dummy",
            "--mq_protein_groups",
            picked_fdr_file,
            "--protein_groups_out",
            "dummy",
            "--fasta",
            fasta_file,
            "--enzyme",
            "trypsinp",
            "--min-length",
            "6",
            "--cleavages",
            "3",
            "--lfq_min_peptide_ratios",
            "1",
            "--gene_level",
        ]
    )

    # obtain map from peptides to proteins to remap the peptide sequences
    parseId = protein_annotation.parse_gene_name_func
    peptideToProteinMaps, _ = quant.get_peptide_to_protein_maps(args, parse_id=parseId)
    peptideToProteinMap = peptideToProteinMaps[0]

    # read in the protein groups from the protein grouping step
    proteinGroupResults = maxquant.parse_mq_protein_groups_file(picked_fdr_file)

    # create a ProteinGroups object that maps protein IDs to their indices in proteinGroupResults
    proteinGroups = ProteinGroups.from_protein_group_results(proteinGroupResults)

    # collect all peptides per protein group per sample
    proteinGroupResults = add_peptides_to_protein_groups(
        df, peptideToProteinMap, proteinGroupResults, proteinGroups
    )
    del df, peptideToProteinMap, proteinGroups

    protein_annotations = protein_annotation.get_protein_annotations_single(
        fasta_file, db="concat", parse_id=parseId
    )

    output_columns: List[columns.ProteinGroupColumns] = [
        columns.ProteinAnnotationsColumns(protein_annotations),
        columns.UniquePeptideCountColumns(),
        columns.LFQIntensityColumns(
            minPeptideRatiosLFQ=args.lfq_min_peptide_ratios,
            stabilizeLargeRatiosLFQ=args.lfq_stabilize_large_ratios,
            numThreads=fdr_num_threads,
        ),
    ]

    for c in output_columns:
        c.append(proteinGroupResults, post_err_prob_cutoff=1.0)

    logger.info(f"Writing picked group results to file")
    proteinGroupResults.write(picked_fdr_file_with_quant)

    logger.info(
        f"Protein group results have been written to: {picked_fdr_file_with_quant}"
    )


def add_peptides_to_protein_groups(
    df,
    peptideToProteinMap,
    proteinGroupResults: ProteinGroupResults,
    proteinGroups: ProteinGroups,
) -> ProteinGroupResults:
    logger.info("Collecting peptides from dataframe per protein group")
    parsedExperiments = set()
    missingPeptidesInFasta = 0

    df.columns = df.columns.str.replace(" ", "_")

    num_tmt_channels = len(
        df.filter(regex="Reporter_intensity_corrected_", axis="columns").columns
    )
    for row in df.itertuples():  # itertuples is ~5x faster than iterrows
        proteins = digest.get_proteins(
            peptideToProteinMap, helpers.clean_peptide(row.Modified_sequence)
        )

        # removes peptides from proteins not present in the fasta file, this often includes peptides from contaminants
        if len(proteins) == 0:
            missingPeptidesInFasta += 1
            continue

        proteins = helpers.remove_decoy_proteins_from_target_peptides(proteins)
        proteinGroupIdxs = proteinGroups.get_protein_group_idxs(proteins)

        if len(proteinGroupIdxs) != 1:
            continue

        proteinGroupIdx = list(proteinGroupIdxs)[0]

        # consider each TMT channel as a separate experiment to activate the MaxLFQ algorithm (pairwise median ratios)
        for tmt_channel in range(1, num_tmt_channels + 1):
            experiment = f"Reporter intensity corrected {tmt_channel} {row.Batch}"
            if experiment not in parsedExperiments:
                parsedExperiments.add(experiment)

            intensity = getattr(row, f"Reporter_intensity_corrected_{tmt_channel}")
            if intensity <= 0.0 or np.isnan(intensity):
                continue

            precursorQuant = PrecursorQuant(
                row.Modified_sequence,
                row.Charge,
                experiment,
                row.Fraction,
                intensity,
                row.PEP,
                None,
                None,
                row.id,
            )
            proteinGroupResults[proteinGroupIdx].precursorQuants.append(precursorQuant)

    if missingPeptidesInFasta > 0:
        logger.info(
            f"Skipped {missingPeptidesInFasta} precursors not present in the fasta file"
        )

    if len(parsedExperiments) > 0:
        proteinGroupResults.experiments = sorted(list(parsedExperiments))

    return proteinGroupResults


def remap_gene_names(df: pd.DataFrame, fasta_file: str) -> pd.DataFrame:
    """
    Re-map gene names based on uniprot identifiers in a fasta file. This is necessary because MaxQuant
    uses their own uniprot->gene mapping file that cannot be changed.
    """
    parseId = protein_annotation.parse_uniprot_id
    proteinAnnotations = protein_annotation.get_protein_annotations_single(
        fasta_file, db="concat", parse_id=parseId
    )

    def get_gene_names(protein_ids: str) -> str:
        """Generates a semicolon separated list of gene names from a semicolon separated list of uniprot protein IDs

        Args:
            protein_ids: semicolon separated list of uniprot protein IDs

        Returns:
            semicolon separated list of gene names
        """
        gene_names = [
            proteinAnnotations.get(
                protein_id,
                protein_annotation.ProteinAnnotation(
                    id=protein_id, fasta_header=protein_id
                ),
            ).gene_name
            for protein_id in protein_ids.split(";")
        ]
        gene_names = {g for g in gene_names if g is not None}
        gene_names_string = ";".join(sorted(gene_names))
        return gene_names_string

    old_gene_names = df["Gene names"]
    df["Gene names"] = df["Proteins"].fillna("").map(get_gene_names)

    # When none of the protein identifiers had a gene name, use the original gene name
    missing_gene_names = df["Gene names"].str.len() == 0
    df.loc[missing_gene_names, "Gene names"] = old_gene_names[missing_gene_names]

    return df


def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-c", "--config", required=True, help="Absolute path to configuration file."
    )
    parser.add_argument(
        "-d",
        "--data_type",
        default="fp",
        help="fp for full proteome and pp for phosphoproteome.",
    )

    args = parser.parse_args(argv)

    configs = config.load(args.config)

    jsonString = json.dumps(configs, indent=4)
    with open(f'{configs["results_folder"]}/configs.json', "w") as jsonFile:
        jsonFile.write(jsonString)

    utils.init_file_logger(configs["results_folder"], "Picked_group_fdr_log.txt")

    sample_annotation_df = sample_annotation.load_sample_annotation(
        configs["sample_annotation"]
    )

    df = pre.load_sample_data(
        configs["results_folder"],
        sample_annotation_df,
        configs["simsi"]["simsi_folder"],
        configs["preprocessing"]["raw_data_location"],
        configs["preprocessing"]["run_simsi"],
        configs["preprocessing"]["run_lfq"],
        configs["preprocessing"]["debug"],
        args.data_type,
    )

    evidence_file_base = "evidence"
    output_file_base = "pickedGeneGroups"
    if args.data_type == "pp":
        df["Modified sequence"] = df["Modified sequence"].str.replace(
            re.compile(r"p([STY])"), lambda pat: f"{pat.group(1)}(Phospho (STY))"
        )
        df = df[df["Modifications"].str.contains("Phospho (STY)", regex=False)]

        evidence_file_base = "evidence_pp"
        output_file_base = "pickedGeneGroups_pp"

    picked_protein_grouping(
        df,
        configs["results_folder"],
        configs["preprocessing"]["picked_fdr"],
        configs["preprocessing"]["fasta_file"],
        configs["preprocessing"]["fdr_num_threads"],
        evidence_file_base,
        output_file_base,
    )


"""
python3 -m topas_pipeline.picked_group -c config_patients.json
"""
if __name__ == "__main__":
    main(sys.argv[1:])
