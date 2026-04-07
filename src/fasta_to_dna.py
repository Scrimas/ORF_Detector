from __future__ import annotations
import re


def parse_fasta_string(fasta_content: str) -> dict[str, str]:
    """
    Parses a FASTA string and extracts the DNA sequences.

    Args:
        fasta_content (str): The FASTA formatted content.

    Returns:
        dict[str, str]: A dictionary mapping sequence IDs to their corresponding DNA sequences.

    Raises:
        ValueError: If a sequence contains non-IUPAC characters.
    """
    sequences: dict[str, list[str]] = {}
    current_id: str = "Unnamed_Sequence"
    iupac_pattern: re.Pattern = re.compile(r"[^A-Z*.-]")

    for line in fasta_content.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            current_id = line[1:].split()[0]
            sequences[current_id] = []
        else:
            if current_id not in sequences:
                sequences[current_id] = []
            cleaned_line: str = line.upper()

            if iupac_pattern.search(cleaned_line):
                invalid_chars: set[str] = set(iupac_pattern.findall(cleaned_line))
                raise ValueError(
                    f"Sequence '{current_id}' contains invalid non-IUPAC characters: {invalid_chars}"
                )
            sequences[current_id].append(cleaned_line)

    result: dict[str, str] = {}
    for seq_id, seq_parts in sequences.items():
        result[seq_id] = "".join(seq_parts)
    return result


def fasta_to_dna(fasta_file: str) -> dict[str, str]:
    """
    Parses a FASTA file and extracts the DNA sequences.

    Args:
        fasta_file (str): Path to the FASTA file.

    Returns:
        dict[str, str]: A dictionary mapping sequence IDs to their corresponding DNA sequences.

    Raises:
        ValueError: If a sequence contains non-IUPAC characters.
    """
    with open(fasta_file, "r", encoding="utf-8") as f:
        fasta_content = f.read()
    return parse_fasta_string(fasta_content)
