def export_orfs_to_txt(orfs_list, output_file_path):
    with open(output_file_path, 'w', encoding='utf-8') as f:
        f.write("=== ANALYSIS REPORT: ORF DETECTION ===\n")
        f.write(f"Total valid ORFs found: {len(orfs_list)}\n")
        f.write("=" * 38 + "\n\n")
        for i, orf in enumerate(orfs_list, start=1):
            seq_id = orf.get("sequence_id", "Unknown")
            strand = orf.get("strand", "Unspecified")
            start = orf.get("start_position", "Unknown")
            end = orf.get("end_position", "Unknown")
            prot_3l = orf.get("protein_3l", "Untranslated")
            prot_1l = orf.get("protein_1l", "Untranslated")
            aa_length = len(prot_1l) if prot_1l != "Untranslated" else 0
            rna = orf.get("rna", "Untranscribed")
            text_block = f"""--- ORF {i} ---
Original Sequence: {seq_id}
Strand: {strand}
Position (+ reference): Nucleotides {start} to {end}
Length: {aa_length} amino acids
RNA Sequence        : {rna}
Sequence (3-letter) : {prot_3l}
Sequence (1-letter) : {prot_1l}

"""
            f.write(text_block)
    print(f"Saved: {output_file_path.name}")