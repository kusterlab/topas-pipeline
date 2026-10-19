AMINO_ACID_CODES = {
    "A": "Alanine",
    "C": "Cysteine",
    "D": "Aspartic Acid",
    "E": "Glutamic Acid",
    "F": "Phenylalanine",
    "G": "Glycine",
    "H": "Histidine",
    "I": "Isoleucine",
    "K": "Lysine",
    "L": "Leucine",
    "M": "Methionine",
    "N": "Asparagine",
    "P": "Proline",
    "Q": "Glutamine",
    "R": "Arginine",
    "S": "Serine",
    "T": "Threonine",
    "V": "Valine",
    "W": "Tryptophan",
    "Y": "Tyrosine",
}
QC_PHOSPHO_SPIKEIN_BACKBONES = {
    "Pool_1": "SLGEHEGGGK",
    "Pool_2": "SSGYDGTSYK",
    "Pool_3": "YASDFGSSAK",
    "Pool_4": "YSYYSSTESK",
    "Pool_5": "YTGAGYDEVK",
    "Pool_6": "SGFSDISIYK",
    "Pool_7": "SDLLSGYYIK",
    "Pool_9": "SIDSFIDSYK",
}
QC_PHOSPHO_SPIKEIN_PEPTIDE_POOLS = {
    f"Pool_{i}": [ "_" + aa + "p" + QC_PHOSPHO_SPIKEIN_BACKBONES[f"Pool_{i}"] + "_" for aa in AMINO_ACID_CODES.keys() ]
    for i in range(1, 10)
    if i != 8
}
QC_PHOSPHO_SPIKEIN_PEPTIDES = [
    peptides
    for pool in QC_PHOSPHO_SPIKEIN_PEPTIDE_POOLS.values()
    for peptides in pool
]