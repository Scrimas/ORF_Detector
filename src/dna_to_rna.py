from __future__ import annotations


def transcribe_dna_to_rna(dna_sequence: str) -> str:
    return dna_sequence.upper().replace("T", "U")
