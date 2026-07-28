MQ_INPUT_COLUMNS = [
    "Modified sequence",
    "Charge",
    "Proteins",
    "Leading proteins",
    "Gene names",
    "Intensity",
    "Raw file",
    "id",
    "Experiment",
    "Type",
    "MS/MS scan number",
    "PEP",
    "Score",
    "Modifications",
    "Reverse",
    "Potential contaminant",
    "Fraction",
]

MQ_OUTPUT_COLUMNS = MQ_INPUT_COLUMNS + [
    "Reporter intensity corrected 1",
    "Batch",
]

DIANN_REPORT_INPUT_COLUMNS = [
    "Modified.Sequence",
    "Precursor.Charge",
    "PEP",
    "Run",
    "Ms1.Normalised",
]

DIANN_REPORT_OUTPUT_COLUMNS = [
    "Modified sequence",
    "Charge",
    "Raw file",
    "id",
    "Intensity",
    "Type",
    "Fraction",
    "PEP",
    "Experiment",
]

DIANN_PSM_INPUT_COLUMNS = [
    "Modified Peptide",
    "Peptide",
    "Charge",
    "Is Decoy",
    "Is Contaminant",
    "Hyperscore",
    "Spectrum",
    "Assigned Modifications",
    "Protein",
    "Mapped Proteins",
    "Gene",
    "Mapped Genes",
]

DIANN_PSM_OUTPUT_COLUMNS = [
    "Modified sequence",
    "Charge",
    "Reverse",
    "Potential contaminant",
    "Score",
    "MS/MS scan number",
    "Modifications",
    "Leading proteins",
    "Proteins",
    "Gene names",
]

IONQUANT_COMBINED_ION_INPUT_COLUMNS = [
    "Modified Sequence",
    "Charge",
    "Protein",
    "Mapped Proteins",
    "Gene",
    "Mapped Genes",
]

IONQUANT_COMBINED_ION_OUTPUT_COLUMNS = [
    "Modified sequence",
    "Charge",
    "Leading proteins",
    "Proteins",
    "Leading gene",
    "Gene names",
    "Intensity",
    "Raw file",
    "id",
    "Type",
    "Experiment",
]

IONQUANT_PSM_INPUT_COLUMNS = DIANN_PSM_INPUT_COLUMNS + ["Probability"]

IONQUANT_PSM_OUTPUT_COLUMNS = [
    "Modified sequence",
    "Charge",
    "Reverse",
    "Potential contaminant",
    "PEP",
    "Score",
    "MS/MS scan number",
    "Modifications",
]


DATA_ACQUISITION_STRATEGY = ["DDA", "DIA"]

QUANT_STRATEGY = ["TMT", "LFQ"]

QUANT_FILE_FORMAT = ["evidence", "ionquant", "diann"]

PSM_TSV_MOD_DICT = {
    "S[167]": "S(Phospho (STY))",
    "T[181]": "T(Phospho (STY))",
    "Y[243]": "Y(Phospho (STY))",
    "M[147]": "M(Oxidation (M))",
    "n[43]": "(Acetyl (Protein N-term))",
}

COMBINED_ION_TSV_MOD_DICT = {
    "S[79.9663]": "S(Phospho (STY))",
    "T[79.9663]": "T(Phospho (STY))",
    "Y[79.9663]": "Y(Phospho (STY))",
    "M[15.9949]": "M(Oxidation (M))",
    "n[42.0106]": "(Acetyl (Protein N-term))",
    "C[57.0215]": "C",
}

COMBINED_ION_TSV_SHORT_MOD_DICT = {
    "S[79.9663]": "S(ph)",
    "T[79.9663]": "T(ph)",
    "Y[79.9663]": "Y(ph)",
    "M[15.9949]": "M(ox)",
    "n[42.0106]": "(ac)",
    "C[57.0215]": "C",
}

DIANN_MOD_DICT = {
    "S(UniMod:21)": "S(Phospho (STY))",
    "T(UniMod:21)": "T(Phospho (STY))",
    "Y(UniMod:21)": "Y(Phospho (STY))",
    "M(UniMod:35)": "M(Oxidation (M))",
    "(UniMod:1)": "(Acetyl (Protein N-term))",
    "C(UniMod:4)": "C",
}

MQ_INPUT_COLUMNS_TMT = MQ_INPUT_COLUMNS + [f"Reporter intensity corrected {i}" for i in range(1, 12)]

MQ_OUTPUT_COLUMNS_TMT = MQ_INPUT_COLUMNS_TMT + ["Batch"]