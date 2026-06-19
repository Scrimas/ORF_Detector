from __future__ import annotations

import re


def parse_fasta_content(fasta_content: str) -> dict[str, str]:
    sequences_accumulator: dict[str, list[str]] = {}
    current_sequence_id: str = "Unnamed_Sequence"
    invalid_character_pattern: re.Pattern = re.compile(r"[^A-Z*.-]")

    for line in fasta_content.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            current_sequence_id = line[1:].split()[0]
            sequences_accumulator[current_sequence_id] = []
        else:
            if current_sequence_id not in sequences_accumulator:
                sequences_accumulator[current_sequence_id] = []
            uppercase_line: str = line.upper()

            if invalid_character_pattern.search(uppercase_line):
                unsupported_characters: set[str] = set(invalid_character_pattern.findall(uppercase_line))
                raise ValueError(
                    f"Sequence '{current_sequence_id}' contains invalid non-IUPAC characters: {unsupported_characters}"
                )
            sequences_accumulator[current_sequence_id].append(uppercase_line)

    parsed_sequences: dict[str, str] = {}
    for sequence_id, sequence_parts in sequences_accumulator.items():
        parsed_sequences[sequence_id] = "".join(sequence_parts)
    return parsed_sequences


def read_fasta_file(fasta_file_path: str) -> dict[str, str]:
    with open(fasta_file_path, "r", encoding="utf-8") as f:
        fasta_content = f.read()
    return parse_fasta_content(fasta_content)
