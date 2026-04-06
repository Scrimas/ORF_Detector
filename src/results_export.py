from __future__ import annotations
from typing import Any

def format_sequence_ncbi(sequence: str, line_length: int = 50, chunk_size: int = 10) -> str:
    """
    Formats a sequence into NCBI-style blocks with position numbers.
    Example:
            1 MALWMRLLPL LALLALWGPD PAAAFVNQHL CGSHLVEALY LVCGERGFFY
           51 TPKTRREAED LQVGQVELGG GPGAGSLQPL ALEGSLQKRG IVEQCCTSIC
    """
    formatted_lines: list[str] = []
    
    # Iterate over the sequence in chunks of `line_length` (default 50)
    for i in range(0, len(sequence), line_length):
        line_seq: str = sequence[i:i + line_length]
        
        # Sub-divide the line into blocks of `chunk_size` (default 10)
        chunks: list[str] = [line_seq[j:j + chunk_size] for j in range(0, len(line_seq), chunk_size)]
        chunked_str: str = " ".join(chunks)
        
        # Right-align the starting position number to 9 characters for the NCBI margin
        formatted_lines.append(f"{i + 1:>9} {chunked_str}")
        
    return "\n".join(formatted_lines)

def export_protein_to_txt(protein_data: dict[str, Any], output_file: str) -> None:
    """
    Exports Protein analysis results to a structured text file.
    """
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write("=" * 60 + "\n")
        f.write("SeqProfiler: Protein Analysis Report\n")
        f.write("=" * 60 + "\n\n")
        
        f.write(f"> {protein_data['sequence_id']}\n")
        f.write("-" * 60 + "\n")
        
        pp: dict[str, Any] = protein_data['prot_props']
        f.write(f"Protein Length: {protein_data['length']} aa  | Mass: {pp['mass_kda']:.2f} kDa  | pI: {pp['pi']:.2f}  | Ext.C: {pp['ext_coeff']} M⁻¹·cm⁻¹\n\n")
        
        f.write("Sequence:\n")
        formatted_seq: str = format_sequence_ncbi(protein_data['protein_1l'])
        f.write(f"{formatted_seq}\n")

def export_orfs_to_txt(orfs: list[dict[str, Any]], output_file: str) -> None:
    """
    Exports ORF analysis results to a structured text file in Classic Bioinformatics format.
    
    Args:
        orfs (List[Dict]): List of ORF data dictionaries.
        output_file (str): Path to the output file.
    """
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write("=" * 60 + "\n")
        f.write("SeqProfiler: ORF Analysis Report\n")
        f.write("=" * 60 + "\n\n")
        
        if not orfs:
            f.write("No ORFs found matching the criteria.\n")
            return

        for orf in orfs:
            # Header
            f.write(f"> {orf['sequence_id']} | Strand: {orf['strand']} | Pos: {orf['start_position']} - {orf['end_position']} bp\n")
            f.write("-" * 60 + "\n")
            
            # Aligned DNA Properties
            dp: dict[str, Any] = orf['dna_props']
            f.write(f"DNA     Length: {dp['length']} bp  | Mass: {dp['mass_da']:,.0f} Da | AT: {dp['at_prop']:.1f}% | GC: {dp['gc_prop']:.1f}% | Tm: {dp['tm']:.1f} °C\n")
            
            # RNA Properties
            rp: dict[str, Any] = orf['rna_props']
            f.write(f"RNA     Mass: {rp['mass_da']:,.0f} Da\n")
            
            # Protein Properties
            pp: dict[str, Any] = orf['prot_props']
            prot_len: int = len(orf['protein_1l'])
            f.write(f"Protein Length: {prot_len} aa  | Mass: {pp['mass_kda']:.2f} kDa  | pI: {pp['pi']:.2f}  | Ext.C: {pp['ext_coeff']} M⁻¹·cm⁻¹\n\n")
            
            # NCBI Chunked Sequence
            f.write("Sequence:\n")
            formatted_seq: str = format_sequence_ncbi(orf['protein_1l'])
            f.write(f"{formatted_seq}\n\n\n\n")