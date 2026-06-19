from __future__ import annotations

import argparse
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any

from dna_to_codon import find_open_reading_frames
from dna_to_protein import translate_sequence
from dna_to_rna import transcribe_dna_to_rna
from fasta_to_dna import parse_fasta_content, read_fasta_file
from ncbi_fetch import fetch_ncbi_sequence
from results_export import export_orf_report, export_protein_report
from sequence_properties import (
    compute_dna_properties,
    compute_protein_properties,
    compute_rna_properties,
)


def display_progress_bar(
    iteration: int,
    total: int,
    prefix: str = "",
    suffix: str = "",
    length: int = 50,
    fill: str = "█",
) -> None:
    percent: str = ("{0:.1f}").format(100 * (iteration / float(total)))
    filled_length: int = int(length * iteration // total)
    bar: str = fill * filled_length + "-" * (length - filled_length)
    sys.stdout.write(f"\r{prefix} |{bar}| {percent}% {suffix}")
    sys.stdout.flush()
    if iteration == total:
        sys.stdout.write("\n")


def process_sequences(
    id_to_sequence_map: dict[str, str],
    dataset_name: str,
    min_orf_length: int,
    allowed_start_codons: set[str],
    output_directory: Path,
    forced_sequence_type: str | None = None,
) -> str:
    try:
        all_orfs: list[dict[str, int | str | Any]] = []
        all_proteins: list[dict[str, Any]] = []
        is_protein_batch = False

        for sequence_id, raw_sequence in id_to_sequence_map.items():
            normalized_sequence = raw_sequence.upper()

            protein_specific_amino_acids = set("EFILPQZJ")
            has_protein_specific_residues = any(char in protein_specific_amino_acids for char in normalized_sequence)

            valid_nucleotides = set("ATGCNURYSWKMBDHV-")
            non_nucleotide_fraction = sum(1 for char in normalized_sequence if char not in valid_nucleotides) / len(normalized_sequence) if len(normalized_sequence) > 0 else 1.0

            if forced_sequence_type == "protein":
                is_protein_sequence = True
            elif forced_sequence_type in ("genomic", "transcript"):
                is_protein_sequence = False
            else:
                is_protein_sequence = has_protein_specific_residues or (non_nucleotide_fraction > 0.1)

            if is_protein_sequence:
                is_protein_batch = True
                protein_properties = compute_protein_properties(raw_sequence)
                protein_record = {
                    "sequence_id": sequence_id,
                    "length": len(raw_sequence.replace("-", "").replace("*", "")),
                    "prot_props": protein_properties,
                    "protein_1l": raw_sequence,
                }
                all_proteins.append(protein_record)
                continue

            sequence_length: int = len(raw_sequence)

            transcript_prefixes = ("NM_", "NR_", "XM_", "XR_")

            if forced_sequence_type == "transcript":
                is_transcript_sequence = True
            elif forced_sequence_type == "genomic":
                is_transcript_sequence = False
            else:
                is_transcript_sequence = sequence_id.startswith(transcript_prefixes)

            forward_orfs: list[dict[str, int | str | Any]] = find_open_reading_frames(
                raw_sequence, min_orf_length=min_orf_length, allowed_start_codons=allowed_start_codons
            )
            for orf in forward_orfs:
                orf.update(
                    {
                        "sequence_id": sequence_id,
                        "strand": "Forward",
                        "rna": transcribe_dna_to_rna(orf["sequence"]),
                    }
                )
                protein_3letter, protein_1letter = translate_sequence(orf["sequence"], is_open_reading_frame=True)
                orf["protein_3l"], orf["protein_1l"] = protein_3letter, protein_1letter
                orf["dna_props"] = compute_dna_properties(orf["sequence"])
                orf["rna_props"] = compute_rna_properties(orf["rna"])
                orf["prot_props"] = compute_protein_properties(protein_1letter)

            if not is_transcript_sequence:
                reverse_complement: str = normalized_sequence.translate(
                    str.maketrans("ATCGU", "TAGCA")
                )[::-1]
                reverse_orfs: list[dict[str, Any]] = find_open_reading_frames(
                    reverse_complement, min_orf_length=min_orf_length, allowed_start_codons=allowed_start_codons
                )
                for orf in reverse_orfs:
                    corrected_start_position: int = sequence_length - orf["start_position"] + 1
                    corrected_end_position: int = sequence_length - orf["end_position"] + 1

                    orf.update(
                        {
                            "sequence_id": sequence_id,
                            "strand": "Reverse",
                            "start_position": corrected_start_position,
                            "end_position": corrected_end_position,
                            "rna": transcribe_dna_to_rna(orf["sequence"]),
                        }
                    )
                    protein_3letter, protein_1letter = translate_sequence(orf["sequence"], is_open_reading_frame=True)
                    orf["protein_3l"], orf["protein_1l"] = protein_3letter, protein_1letter
                    orf["dna_props"] = compute_dna_properties(orf["sequence"])
                    orf["rna_props"] = compute_rna_properties(orf["rna"])
                    orf["prot_props"] = compute_protein_properties(protein_1letter)

                all_orfs.extend(forward_orfs + reverse_orfs)
            else:
                all_orfs.extend(forward_orfs)

        if is_protein_batch:
            report_file_path: Path = output_directory / f"results_{dataset_name}.txt"
            export_protein_report(all_proteins, str(report_file_path))
            return f"Successfully processed protein batch {dataset_name}"

        all_orfs = sorted(
            all_orfs, key=lambda x: (x["sequence_id"], x["start_position"])
        )
        report_file_path: Path = output_directory / f"results_{dataset_name}.txt"
        export_orf_report(all_orfs, str(report_file_path))

        return f"Successfully processed {dataset_name}"
    except Exception as e:
        return f"Error processing {dataset_name}: {str(e)}"


def process_single_file(
    input_file_path: Path,
    min_orf_length: int,
    allowed_start_codons: set[str],
    output_directory: Path,
    forced_sequence_type: str | None = None,
) -> str:
    try:
        id_to_sequence_map: dict[str, str] = read_fasta_file(str(input_file_path))
        return process_sequences(
            id_to_sequence_map,
            input_file_path.stem,
            min_orf_length,
            allowed_start_codons,
            output_directory,
            forced_sequence_type,
        )
    except Exception as e:
        return f"Error reading {input_file_path.name}: {str(e)}"


def process_ncbi_id(
    accession_number: str,
    min_orf_length: int,
    allowed_start_codons: set[str],
    output_directory: Path,
    forced_sequence_type: str | None = None,
) -> str:
    try:
        fetched_fasta = fetch_ncbi_sequence(accession_number)
        id_to_sequence_map = parse_fasta_content(fetched_fasta)
        return process_sequences(
            id_to_sequence_map,
            accession_number,
            min_orf_length,
            allowed_start_codons,
            output_directory,
            forced_sequence_type,
        )
    except Exception as e:
        return f"Error processing NCBI ID {accession_number}: {str(e)}"


def main() -> None:
    project_root: Path = Path(__file__).resolve().parent.parent

    arg_parser: argparse.ArgumentParser = argparse.ArgumentParser(
        description="SeqProfiler: High-performance DNA analysis tool."
    )
    arg_parser.add_argument(
        "--min-length",
        type=int,
        default=50,
        help="Minimum ORF size (in amino acids) [default: 50]",
    )
    arg_parser.add_argument(
        "--input",
        type=str,
        default=None,
        help="Path to input directory [default: data/]",
    )
    arg_parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Path to output directory [default: results/]",
    )
    arg_parser.add_argument(
        "--workers",
        type=int,
        default=None,
        help="Number of parallel workers [default: CPU count]",
    )
    arg_parser.add_argument(
        "--start-codons",
        type=str,
        default="ATG",
        help="Comma-separated list of alternative start codons (e.g., ATG,CTG,GTG) [default: ATG]",
    )
    arg_parser.add_argument(
        "--ncbi",
        type=str,
        default=None,
        help="Comma-separated list of NCBI accession IDs to fetch and analyze",
    )
    arg_parser.add_argument(
        "--seq-type",
        type=str,
        choices=["auto", "genomic", "transcript", "protein"],
        default="auto",
        help="Force sequence type. 'genomic' (dsDNA, search both strands), 'transcript' (mRNA/cDNA, search forward only), or 'protein'. [default: auto]",
    )
    parsed_args: argparse.Namespace = arg_parser.parse_args()

    allowed_start_codons: set[str] = {
        codon.strip().upper() for codon in parsed_args.start_codons.split(",") if codon.strip()
    }

    forced_sequence_type = parsed_args.seq_type if parsed_args.seq_type != "auto" else None

    print("\n" + "=" * 50)
    print(" " * 19 + "SeqProfiler")
    print("=" * 50 + "\n")

    is_using_default_data = parsed_args.input is None and parsed_args.ncbi is None
    input_directory: Path | None = None
    if parsed_args.input:
        input_directory = Path(parsed_args.input).resolve()
    elif is_using_default_data:
        input_directory = project_root / "data"

    output_directory: Path = (
        Path(parsed_args.output).resolve() if parsed_args.output else project_root / "results"
    )
    output_directory.mkdir(parents=True, exist_ok=True)

    analysis_tasks = []

    if input_directory and input_directory.exists():
        input_fasta_files: list[Path] = list(input_directory.glob("*.fasta"))
        for f in input_fasta_files:
            analysis_tasks.append(("file", f))
    elif input_directory and not input_directory.exists():
        print(f"Error: Input directory '{input_directory}' does not exist.")
        sys.exit(1)

    if parsed_args.ncbi:
        accession_numbers = [nid.strip() for nid in parsed_args.ncbi.split(",") if nid.strip()]
        for nid in accession_numbers:
            analysis_tasks.append(("ncbi", nid))

    if not analysis_tasks:
        print("No input files or NCBI IDs provided. Exiting.")
        return

    source_description = "None"
    if parsed_args.input and parsed_args.ncbi:
        source_description = f"Files ({input_directory}) + NCBI"
    elif parsed_args.input:
        source_description = f"Files ({input_directory})"
    elif parsed_args.ncbi:
        source_description = "NCBI"
    elif is_using_default_data:
        source_description = f"Default Files ({input_directory})"

    print(f"[*] Configuration:")
    print(f"    - Input:        {source_description}")
    print(f"    - Output:       {output_directory}")
    print(f"    - Seq Type:     {parsed_args.seq_type}")
    print(f"    - Start Codons: {', '.join(allowed_start_codons)}")
    print(f"    - Min AA:       {parsed_args.min_length}")
    print(f"    - Total tasks:  {len(analysis_tasks)}\n")

    print(f"[*] Starting Analysis...")
    display_progress_bar(0, len(analysis_tasks), prefix="Progress:", suffix="Complete", length=40)

    results: list[str] = []
    with ProcessPoolExecutor(max_workers=parsed_args.workers) as executor:
        futures = []
        for task_type, task_val in analysis_tasks:
            if task_type == "file":
                futures.append(
                    executor.submit(
                        process_single_file,
                        task_val,
                        parsed_args.min_length,
                        allowed_start_codons,
                        output_directory,
                        forced_sequence_type,
                    )
                )
            else:
                futures.append(
                    executor.submit(
                        process_ncbi_id,
                        task_val,
                        parsed_args.min_length,
                        allowed_start_codons,
                        output_directory,
                        forced_sequence_type,
                    )
                )

        completed_tasks: int = 0
        for future in as_completed(futures):
            results.append(future.result())
            completed_tasks += 1
            display_progress_bar(
                completed_tasks, len(analysis_tasks), prefix="Progress:", suffix="Complete", length=40
            )

    print("\n[*] Summary:")
    for res in results:
        print(f"    - {res}")

    print("\nAnalysis finished successfully.")


if __name__ == "__main__":
    main()
