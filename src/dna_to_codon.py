from __future__ import annotations


def find_open_reading_frames(
    dna_sequence: str,
    min_orf_length: int = 50,
    allowed_start_codons: set[str] | None = None,
) -> list[dict[str, int | str]]:
    if allowed_start_codons is None:
        allowed_start_codons = {"ATG"}

    allowed_start_codons = {codon.replace("U", "T") for codon in allowed_start_codons}

    found_open_reading_frames: list[dict[str, int | str]] = []
    normalized_sequence: str = dna_sequence.upper().replace("U", "T")
    sequence_length: int = len(normalized_sequence)
    termination_codons: set[str] = {"TAA", "TAG", "TGA"}
    tracked_stop_positions: set[int] = set()

    for reading_frame in range(3):
        for i in range(reading_frame, sequence_length - 2, 3):
            codon: str = normalized_sequence[i : i + 3]

            if codon in allowed_start_codons:
                for j in range(i, sequence_length - 2, 3):
                    candidate_codon: str = normalized_sequence[j : j + 3]
                    if candidate_codon in termination_codons:
                        end_position: int = j + 3

                        if end_position not in tracked_stop_positions:
                            orf_sequence: str = dna_sequence[i:end_position]
                            if len(orf_sequence) >= (min_orf_length * 3):
                                found_open_reading_frames.append(
                                    {
                                        "start_position": i + 1,
                                        "end_position": end_position,
                                        "sequence": orf_sequence,
                                    }
                                )
                                tracked_stop_positions.add(end_position)
                        break
    return found_open_reading_frames
