from __future__ import annotations

CODON_TO_AMINO_ACID: dict[str, str] = {
    "AAA": "Lys",
    "AAT": "Asn",
    "AAG": "Lys",
    "AAC": "Asn",
    "ATA": "Ile",
    "ATT": "Ile",
    "ATG": "Met",
    "ATC": "Ile",
    "AGA": "Arg",
    "AGT": "Ser",
    "AGG": "Arg",
    "AGC": "Ser",
    "ACA": "Thr",
    "ACT": "Thr",
    "ACG": "Thr",
    "ACC": "Thr",
    "TAA": "Ter",
    "TAT": "Tyr",
    "TAG": "Ter",
    "TAC": "Tyr",
    "TTA": "Leu",
    "TTT": "Phe",
    "TTG": "Leu",
    "TTC": "Phe",
    "TGA": "Ter",
    "TGT": "Cys",
    "TGG": "Trp",
    "TGC": "Cys",
    "TCA": "Ser",
    "TCT": "Ser",
    "TCG": "Ser",
    "TCC": "Ser",
    "GAA": "Glu",
    "GAT": "Asp",
    "GAG": "Glu",
    "GAC": "Asp",
    "GTA": "Val",
    "GTT": "Val",
    "GTG": "Val",
    "GTC": "Val",
    "GGA": "Gly",
    "GGT": "Gly",
    "GGG": "Gly",
    "GGC": "Gly",
    "GCA": "Ala",
    "GCT": "Ala",
    "GCG": "Ala",
    "GCC": "Ala",
    "CAA": "Gln",
    "CAT": "His",
    "CAG": "Gln",
    "CAC": "His",
    "CTA": "Leu",
    "CTT": "Leu",
    "CTG": "Leu",
    "CTC": "Leu",
    "CGA": "Arg",
    "CGT": "Arg",
    "CGG": "Arg",
    "CGC": "Arg",
    "CCA": "Pro",
    "CCT": "Pro",
    "CCG": "Pro",
    "CCC": "Pro",
}

AMINO_ACID_3_TO_1: dict[str, str] = {
    "Lys": "K",
    "Asn": "N",
    "Ile": "I",
    "Met": "M",
    "Arg": "R",
    "Ser": "S",
    "Thr": "T",
    "Ter": "*",
    "Tyr": "Y",
    "Leu": "L",
    "Phe": "F",
    "Cys": "C",
    "Trp": "W",
    "Glu": "E",
    "Asp": "D",
    "Val": "V",
    "Gly": "G",
    "Ala": "A",
    "Gln": "Q",
    "His": "H",
    "Pro": "P",
    "X": "X",
}


def translate_codon(codon: str) -> str:
    return CODON_TO_AMINO_ACID.get(codon.upper().replace("U", "T"), "X")


def translate_sequence(
    nucleotide_sequence: str, is_open_reading_frame: bool = False
) -> tuple[str, str]:
    protein_3letter_list: list[str] = []
    protein_1letter_list: list[str] = []

    for i in range(0, len(nucleotide_sequence) - 2, 3):
        codon: str = nucleotide_sequence[i : i + 3]

        if i == 0 and is_open_reading_frame:
            amino_acid_3letter = "Met"
        else:
            amino_acid_3letter = translate_codon(codon)

        if amino_acid_3letter == "Ter":
            break

        protein_3letter_list.append(amino_acid_3letter)
        protein_1letter_list.append(AMINO_ACID_3_TO_1.get(amino_acid_3letter, "?"))

    return "-".join(protein_3letter_list), "".join(protein_1letter_list)
