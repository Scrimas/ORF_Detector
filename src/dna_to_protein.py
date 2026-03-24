from __future__ import annotations

GENETIC_CODE: dict[str, str] = {
    "AAA": "Lys", "AAT": "Asn", "AAG": "Lys", "AAC": "Asn",
    "ATA": "Ile", "ATT": "Ile", "ATG": "Met", "ATC": "Ile",
    "AGA": "Arg", "AGT": "Ser", "AGG": "Arg", "AGC": "Ser",
    "ACA": "Thr", "ACT": "Thr", "ACG": "Thr", "ACC": "Thr",
    "TAA": "Ter", "TAT": "Tyr", "TAG": "Ter", "TAC": "Tyr",
    "TTA": "Leu", "TTT": "Phe", "TTG": "Leu", "TTC": "Phe",
    "TGA": "Ter", "TGT": "Cys", "TGG": "Trp", "TGC": "Cys",
    "TCA": "Ser", "TCT": "Ser", "TCG": "Ser", "TCC": "Ser",
    "GAA": "Glu", "GAT": "Asp", "GAG": "Glu", "GAC": "Asp",
    "GTA": "Val", "GTT": "Val", "GTG": "Val", "GTC": "Val",
    "GGA": "Gly", "GGT": "Gly", "GGG": "Gly", "GGC": "Gly",
    "GCA": "Ala", "GCT": "Ala", "GCG": "Ala", "GCC": "Ala",
    "CAA": "Gln", "CAT": "His", "CAG": "Gln", "CAC": "His",
    "CTA": "Leu", "CTT": "Leu", "CTG": "Leu", "CTC": "Leu",
    "CGA": "Arg", "CGT": "Arg", "CGG": "Arg", "CGC": "Arg",
    "CCA": "Pro", "CCT": "Pro", "CCG": "Pro", "CCC": "Pro"
}

AA_3_TO_1: dict[str, str] = {
    "Lys": "K", "Asn": "N", "Ile": "I", "Met": "M", "Arg": "R", "Ser": "S",
    "Thr": "T", "Ter": "*", "Tyr": "Y", "Leu": "L", "Phe": "F", "Cys": "C",
    "Trp": "W", "Glu": "E", "Asp": "D", "Val": "V", "Gly": "G", "Ala": "A",
    "Gln": "Q", "His": "H", "Pro": "P", "X": "X"
}

def dna_to_protein(dna_codon: str) -> str:
    """
    Translates a single DNA codon into its 3-letter amino acid representation.

    Args:
        dna_codon (str): A 3-character DNA codon.

    Returns:
        str: The 3-letter amino acid code (e.g., "Met", "Ser") or "X" if unknown.
    """
    return GENETIC_CODE.get(dna_codon.upper(), "X")

def translate_sequence(dna_sequence: str) -> tuple[str, str]:
    """
    Translates a DNA sequence into protein sequences (both 3-letter and 1-letter formats).

    Args:
        dna_sequence (str): The DNA sequence to translate.

    Returns:
        Tuple[str, str]: A tuple containing (3-letter protein sequence, 1-letter protein sequence).
    """
    protein_3l: list[str] = []
    protein_1l: list[str] = []
    
    for i in range(0, len(dna_sequence) - 2, 3):
        codon: str = dna_sequence[i:i+3]
        aa_3: str = dna_to_protein(codon)
        
        if aa_3 == "Ter":
            break
            
        protein_3l.append(aa_3)
        protein_1l.append(AA_3_TO_1.get(aa_3, "?")) 
        
    return "-".join(protein_3l), "".join(protein_1l)
