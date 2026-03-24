from __future__ import annotations

def dna_to_rna(dna_sequence: str) -> str:
    """
    Transcribes a DNA sequence into an mRNA sequence by replacing Thymine (T) with Uracil (U).

    Args:
        dna_sequence (str): The DNA sequence to transcribe.

    Returns:
        str: The transcribed mRNA sequence.
    """
    return dna_sequence.upper().replace("T", "U")
