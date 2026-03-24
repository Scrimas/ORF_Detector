import argparse
import sys
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import List, Dict, Any, Set

from dna_to_rna import dna_to_rna
from fasta_to_dna import fasta_to_dna
from dna_to_codon import get_orfs
from dna_to_protein import translate_sequence
from results_export import export_orfs_to_txt
from sequence_properties import (
    calculate_dna_properties, 
    calculate_rna_properties, 
    calculate_protein_properties
)

def print_progress_bar(iteration: int, total: int, prefix: str = '', suffix: str = '', length: int = 50, fill: str = '█'):
    """
    Call in a loop to create terminal progress bar
    """
    percent = ("{0:.1f}").format(100 * (iteration / float(total)))
    filled_length = int(length * iteration // total)
    bar = fill * filled_length + '-' * (length - filled_length)
    sys.stdout.write(f'\r{prefix} |{bar}| {percent}% {suffix}')
    sys.stdout.flush()
    if iteration == total:
        sys.stdout.write('\n')

def process_single_file(file_path: Path, min_length: int, start_codons: Set[str], results_dir: Path) -> str:
    """
    Processes a single FASTA file: identifies ORFs, calculates properties, and exports results.
    Designed to be run in parallel.
    """
    try:
        sequences_dict = fasta_to_dna(str(file_path))
        all_file_orfs: List[Dict[str, Any]] = []
        
        for seq_id, main_sequence in sequences_dict.items():
            seq_len = len(main_sequence)
            
            # Forward Strand Processing
            positive_orfs = get_orfs(main_sequence, min_length_aa=min_length, start_codons=start_codons)
            for orf in positive_orfs:
                orf.update({
                    "sequence_id": seq_id,
                    "strand": "Forward",
                    "rna": dna_to_rna(orf["sequence"])
                })
                p3, p1 = translate_sequence(orf["sequence"])
                orf["protein_3l"], orf["protein_1l"] = p3, p1
                orf["dna_props"] = calculate_dna_properties(orf["sequence"])
                orf["rna_props"] = calculate_rna_properties(orf["rna"])
                orf["prot_props"] = calculate_protein_properties(p1)
                
            # Reverse Strand Processing
            reverse_complement = main_sequence.translate(str.maketrans("ATCG", "TAGC"))[::-1]
            negative_orfs = get_orfs(reverse_complement, min_length_aa=min_length, start_codons=start_codons)
            for orf in negative_orfs:

                true_start = seq_len - orf["end_position"] + 1
                true_end = seq_len - orf["start_position"] + 1
                
                orf.update({
                    "sequence_id": seq_id,
                    "strand": "Reverse",
                    "start_position": true_start,
                    "end_position": true_end,
                    "rna": dna_to_rna(orf["sequence"])
                })
                p3, p1 = translate_sequence(orf["sequence"])
                orf["protein_3l"], orf["protein_1l"] = p3, p1
                orf["dna_props"] = calculate_dna_properties(orf["sequence"])
                orf["rna_props"] = calculate_rna_properties(orf["rna"])
                orf["prot_props"] = calculate_protein_properties(p1)
                
            all_file_orfs.extend(positive_orfs + negative_orfs)
            
        all_file_orfs = sorted(all_file_orfs, key=lambda x: (x["sequence_id"], x["start_position"]))
        output_path = results_dir / f"results_{file_path.stem}.txt"
        export_orfs_to_txt(all_file_orfs, str(output_path))
        
        return f"Successfully processed {file_path.name}"
    except Exception as e:
        return f"Error processing {file_path.name}: {str(e)}"

def main():
    base_path = Path(__file__).resolve().parent.parent
    
    parser = argparse.ArgumentParser(description="SeqProfiler: High-performance DNA analysis tool.")
    parser.add_argument("--min-length", type=int, default=50, help="Minimum ORF size (in amino acids) [default: 50]")
    parser.add_argument("--input", type=str, default=None, help="Path to input directory [default: data/]")
    parser.add_argument("--output", type=str, default=None, help="Path to output directory [default: results/]")
    parser.add_argument("--workers", type=int, default=None, help="Number of parallel workers [default: CPU count]")
    parser.add_argument("--start-codons", type=str, default="ATG", help="Comma-separated list of alternative start codons (e.g., ATG,CTG,GTG) [default: ATG]")
    args = parser.parse_args()

    start_codons_set = {codon.strip().upper() for codon in args.start_codons.split(",") if codon.strip()}

    print("\n" + "="*50)
    print(" "*19 + "SeqProfiler")
    print("="*50 + "\n")

    data_dir = Path(args.input).resolve() if args.input else base_path / "data"
    results_dir = Path(args.output).resolve() if args.output else base_path / "results"
    
    if not data_dir.exists():
        print(f"Error: Input directory '{data_dir}' does not exist.")
        sys.exit(1)
    
    results_dir.mkdir(parents=True, exist_ok=True)
    fasta_files = list(data_dir.glob("*.fasta"))
    
    if not fasta_files:
        print("No .fasta files found. Exiting.")
        return

    print(f"[*] Configuration:")
    print(f"    - Input:        {data_dir}")
    print(f"    - Output:       {results_dir}")
    print(f"    - Start Codons: {', '.join(start_codons_set)}")
    print(f"    - Min AA:       {args.min_length}")
    print(f"    - Files:        {len(fasta_files)}\n")

    print(f"[*] Starting Analysis...")
    print_progress_bar(0, len(fasta_files), prefix='Progress:', suffix='Complete', length=40)
    
    results = []
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        future_to_file = {executor.submit(process_single_file, f, args.min_length, start_codons_set, results_dir): f for f in fasta_files}
        
        completed = 0
        for future in as_completed(future_to_file):
            results.append(future.result())
            completed += 1
            print_progress_bar(completed, len(fasta_files), prefix='Progress:', suffix='Complete', length=40)

    print("\n[*] Summary:")
    for res in results:
        print(f"    - {res}")
    
    print("\nAnalysis finished successfully.")

if __name__ == "__main__":
    main()