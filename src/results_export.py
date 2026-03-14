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
            dna_props = orf.get("dna_props", {})
            rna_props = orf.get("rna_props", {})
            prot_props = orf.get("prot_props", {})
            prot_3l = orf.get("protein_3l", "Untranslated")
            prot_1l = orf.get("protein_1l", "Untranslated")
            rna = orf.get("rna", "Untranscribed")
            text_block = f"""--- ORF {i} ---
Original Sequence: {seq_id}
Strand: {strand}
Position (Forward reference): Nucleotides {start} to {end}

[ DNA Properties ]
- Length: {dna_props.get("length", 0)} bp
- Mass: {dna_props.get("mass_da", 0):.2f} Da
- GC Content: {dna_props.get("gc_prop", 0):.2f} %
- AT Content: {dna_props.get("at_prop", 0):.2f} %
- Melting Point (Tm): {dna_props.get("tm", 0):.2f} °C

[ RNA Properties ]
- Length: {rna_props.get("length", 0)} nt
- Mass: {rna_props.get("mass_da", 0):.2f} Da
RNA Sequence:
{rna}

[ Protein Properties ]
- Length: {len(prot_1l) if prot_1l != "Untranslated" else 0} amino acids
- Molecular Weight: {prot_props.get("mass_kda", 0):.2f} kDa
- Theoretical Isoelectric Point (pI): {prot_props.get("pi", 0):.2f}
- Extinction Coefficient (280nm): {prot_props.get("ext_coeff", 0)} M⁻¹ cm⁻¹

Sequence (3-letter):
{prot_3l}
Sequence (1-letter):
{prot_1l}

"""
            f.write(text_block)
    print(f"Saved: {output_file_path.name}")