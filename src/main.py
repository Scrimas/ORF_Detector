from __future__ import annotations
import argparse
import sys
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Any

from dna_to_rna import dna_to_rna
from fasta_to_dna import fasta_to_dna, parse_fasta_string
from ncbi_fetch import fetch_ncbi_sequence
from dna_to_codon import get_orfs
from dna_to_protein import translate_sequence
from results_export import export_orfs_to_txt, export_protein_to_txt
from sequence_properties import (
    calculate_dna_properties,
    calculate_rna_properties,
    calculate_protein_properties,
)


def print_progress_bar(
    iteration: int,
    total: int,
    prefix: str = "",
    suffix: str = "",
    length: int = 50,
    fill: str = "█",
) -> None:
    """
    Call in a loop to create terminal progress bar
    """
    percent: str = ("{0:.1f}").format(100 * (iteration / float(total)))
    filled_length: int = int(length * iteration // total)
    bar: str = fill * filled_length + "-" * (length - filled_length)
    sys.stdout.write(f"\r{prefix} |{bar}| {percent}% {suffix}")
    sys.stdout.flush()
    if iteration == total:
        sys.stdout.write("\n")


def process_sequences(
    sequences_dict: dict[str, str],
    source_name: str,
    min_length: int,
    start_codons: set[str],
    results_dir: Path,
) -> str:
    """
    Processes a dictionary of sequences: identifies ORFs or calculates protein properties.
    """
    try:
        all_file_orfs: list[dict[str, int | str | Any]] = []
        is_protein_batch = False

        for seq_id, main_sequence in sequences_dict.items():
            nucleotides = set("ATGCNU")
            is_protein = any(char not in nucleotides for char in main_sequence.upper())

            if is_protein:
                is_protein_batch = True
                prot_props = calculate_protein_properties(main_sequence)
                protein_data = {
                    "sequence_id": seq_id,
                    "length": len(main_sequence),
                    "prot_props": prot_props,
                    "protein_1l": main_sequence,
                }
                output_path: Path = results_dir / f"results_{source_name}.txt"
                export_protein_to_txt(protein_data, str(output_path))
                continue

            seq_len: int = len(main_sequence)

            # Forward Strand Processing
            positive_orfs: list[dict[str, int | str | Any]] = get_orfs(
                main_sequence, min_length_aa=min_length, start_codons=start_codons
            )
            for orf in positive_orfs:
                orf.update(
                    {
                        "sequence_id": seq_id,
                        "strand": "Forward",
                        "rna": dna_to_rna(orf["sequence"]),
                    }
                )
                p3, p1 = translate_sequence(orf["sequence"])
                orf["protein_3l"], orf["protein_1l"] = p3, p1
                orf["dna_props"] = calculate_dna_properties(orf["sequence"])
                orf["rna_props"] = calculate_rna_properties(orf["rna"])
                orf["prot_props"] = calculate_protein_properties(p1)

            # Reverse Strand Processing
            reverse_complement: str = main_sequence.translate(
                str.maketrans("ATCG", "TAGC")
            )[::-1]
            negative_orfs: list[dict[str, Any]] = get_orfs(
                reverse_complement, min_length_aa=min_length, start_codons=start_codons
            )
            for orf in negative_orfs:

                true_start: int = seq_len - orf["end_position"] + 1
                true_end: int = seq_len - orf["start_position"] + 1

                orf.update(
                    {
                        "sequence_id": seq_id,
                        "strand": "Reverse",
                        "start_position": true_start,
                        "end_position": true_end,
                        "rna": dna_to_rna(orf["sequence"]),
                    }
                )
                p3, p1 = translate_sequence(orf["sequence"])
                orf["protein_3l"], orf["protein_1l"] = p3, p1
                orf["dna_props"] = calculate_dna_properties(orf["sequence"])
                orf["rna_props"] = calculate_rna_properties(orf["rna"])
                orf["prot_props"] = calculate_protein_properties(p1)

            all_file_orfs.extend(positive_orfs + negative_orfs)

        if is_protein_batch:
            return f"Successfully processed protein {source_name}"

        all_file_orfs = sorted(
            all_file_orfs, key=lambda x: (x["sequence_id"], x["start_position"])
        )
        output_path: Path = results_dir / f"results_{source_name}.txt"
        export_orfs_to_txt(all_file_orfs, str(output_path))

        return f"Successfully processed {source_name}"
    except Exception as e:
        return f"Error processing {source_name}: {str(e)}"


def process_single_file(
    file_path: Path, min_length: int, start_codons: set[str], results_dir: Path
) -> str:
    """
    Processes a single FASTA file.
    """
    try:
        sequences_dict: dict[str, str] = fasta_to_dna(str(file_path))
        return process_sequences(
            sequences_dict, file_path.stem, min_length, start_codons, results_dir
        )
    except Exception as e:
        return f"Error reading {file_path.name}: {str(e)}"


def process_ncbi_id(
    ncbi_id: str, min_length: int, start_codons: set[str], results_dir: Path
) -> str:
    """
    Fetches and processes a single NCBI ID.
    """
    try:
        fasta_content = fetch_ncbi_sequence(ncbi_id)
        sequences_dict = parse_fasta_string(fasta_content)
        return process_sequences(
            sequences_dict, ncbi_id, min_length, start_codons, results_dir
        )
    except Exception as e:
        return f"Error processing NCBI ID {ncbi_id}: {str(e)}"


def main() -> None:
    base_path: Path = Path(__file__).resolve().parent.parent

    parser: argparse.ArgumentParser = argparse.ArgumentParser(
        description="SeqProfiler: High-performance DNA analysis tool."
    )
    parser.add_argument(
        "--min-length",
        type=int,
        default=50,
        help="Minimum ORF size (in amino acids) [default: 50]",
    )
    parser.add_argument(
        "--input",
        type=str,
        default=None,
        help="Path to input directory [default: data/]",
    )
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Path to output directory [default: results/]",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=None,
        help="Number of parallel workers [default: CPU count]",
    )
    parser.add_argument(
        "--start-codons",
        type=str,
        default="ATG",
        help="Comma-separated list of alternative start codons (e.g., ATG,CTG,GTG) [default: ATG]",
    )
    parser.add_argument(
        "--ncbi",
        type=str,
        default=None,
        help="Comma-separated list of NCBI accession IDs to fetch and analyze",
    )
    args: argparse.Namespace = parser.parse_args()

    start_codons_set: set[str] = {
        codon.strip().upper() for codon in args.start_codons.split(",") if codon.strip()
    }

    print("\n" + "=" * 50)
    print(" " * 19 + "SeqProfiler")
    print("=" * 50 + "\n")

    use_default_data = args.input is None and args.ncbi is None
    data_dir: Path | None = None
    if args.input:
        data_dir = Path(args.input).resolve()
    elif use_default_data:
        data_dir = base_path / "data"

    results_dir: Path = (
        Path(args.output).resolve() if args.output else base_path / "results"
    )
    results_dir.mkdir(parents=True, exist_ok=True)

    tasks = []

    if data_dir and data_dir.exists():
        fasta_files: list[Path] = list(data_dir.glob("*.fasta"))
        for f in fasta_files:
            tasks.append(("file", f))
    elif data_dir and not data_dir.exists():
        print(f"Error: Input directory '{data_dir}' does not exist.")
        sys.exit(1)

    if args.ncbi:
        ncbi_ids = [nid.strip() for nid in args.ncbi.split(",") if nid.strip()]
        for nid in ncbi_ids:
            tasks.append(("ncbi", nid))

    if not tasks:
        print("No input files or NCBI IDs provided. Exiting.")
        return

    input_source = "None"
    if args.input and args.ncbi:
        input_source = f"Files ({data_dir}) + NCBI"
    elif args.input:
        input_source = f"Files ({data_dir})"
    elif args.ncbi:
        input_source = "NCBI"
    elif use_default_data:
        input_source = f"Default Files ({data_dir})"

    print(f"[*] Configuration:")
    print(f"    - Input:        {input_source}")
    print(f"    - Output:       {results_dir}")
    print(f"    - Start Codons: {', '.join(start_codons_set)}")
    print(f"    - Min AA:       {args.min_length}")
    print(f"    - Total tasks:  {len(tasks)}\n")

    print(f"[*] Starting Analysis...")
    print_progress_bar(0, len(tasks), prefix="Progress:", suffix="Complete", length=40)

    results: list[str] = []
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = []
        for task_type, task_val in tasks:
            if task_type == "file":
                futures.append(
                    executor.submit(
                        process_single_file,
                        task_val,
                        args.min_length,
                        start_codons_set,
                        results_dir,
                    )
                )
            else:
                futures.append(
                    executor.submit(
                        process_ncbi_id,
                        task_val,
                        args.min_length,
                        start_codons_set,
                        results_dir,
                    )
                )

        completed: int = 0
        for future in as_completed(futures):
            results.append(future.result())
            completed += 1
            print_progress_bar(
                completed, len(tasks), prefix="Progress:", suffix="Complete", length=40
            )

    print("\n[*] Summary:")
    for res in results:
        print(f"    - {res}")

    print("\nAnalysis finished successfully.")


if __name__ == "__main__":
    main()
