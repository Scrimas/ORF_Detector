from __future__ import annotations


def get_orfs(
    dna_sequence: str, min_length_aa: int = 50, start_codons: set[str] | None = None
) -> list[dict[str, int | str]]:
    """
    Identifies Open Reading Frames (ORFs) in a DNA sequence.

    Optimized to find non-redundant ORFs (longest frame per stop codon) and uses
    frame-based scanning for better efficiency. Supports alternative start codons.

    Args:
        dna_sequence (str): The DNA sequence to scan.
        min_length_aa (int): Minimum length of the ORF in amino acids.
        start_codons (Set[str], optional): Allowed start codons. Defaults to {"ATG"}.

    Returns:
        List[Dict]: A list of dictionaries containing start, end, and sequence of each ORF.
    """
    if start_codons is None:
        start_codons = {"ATG"}

    found_orfs: list[dict[str, int | str]] = []
    seq_len: int = len(dna_sequence)
    stop_codons: set[str] = {"TAA", "TAG", "TGA"}
    used_stops: set[int] = set()

    for frame in range(3):
        for i in range(frame, seq_len - 2, 3):
            current_codon: str = dna_sequence[i : i + 3]

            if current_codon in start_codons:
                for j in range(i, seq_len - 2, 3):
                    reading_codon: str = dna_sequence[j : j + 3]
                    if reading_codon in stop_codons:
                        end_pos: int = j + 3

                        # Filtering heuristic: only keep the longest ORF for a given stop
                        # This should prevents reporting nested, smaller ORFs within the same frame
                        if end_pos not in used_stops:
                            orf_sequence: str = dna_sequence[i:end_pos]
                            if len(orf_sequence) >= (min_length_aa * 3):
                                found_orfs.append(
                                    {
                                        "start_position": i + 1,
                                        "end_position": end_pos,
                                        "sequence": orf_sequence,
                                    }
                                )
                                used_stops.add(end_pos)
                        break
    return found_orfs
