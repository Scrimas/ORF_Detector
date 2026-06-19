from __future__ import annotations

from typing import Any


def format_sequence_ncbi_style(
    sequence_str: str, line_width: int = 50, block_size: int = 10
) -> str:
    output_lines: list[str] = []

    for i in range(0, len(sequence_str), line_width):
        line_content: str = sequence_str[i : i + line_width]

        line_blocks: list[str] = [
            line_content[j : j + block_size] for j in range(0, len(line_content), block_size)
        ]
        formatted_line_content: str = " ".join(line_blocks)

        output_lines.append(f"{i + 1:>9} {formatted_line_content}")

    return "\n".join(output_lines)


def export_protein_report(
    protein_records: list[dict[str, Any]], report_file_path: str
) -> None:
    with open(report_file_path, "w", encoding="utf-8") as f:
        f.write("=" * 60 + "\n")
        f.write("SeqProfiler: Protein Analysis Report\n")
        f.write("=" * 60 + "\n\n")

        for protein_record in protein_records:
            f.write(f"> {protein_record['sequence_id']}\n")
            f.write("-" * 60 + "\n")

            protein_properties: dict[str, Any] = protein_record["prot_props"]
            f.write(
                f"Protein Length: {protein_record['length']} aa  | Mass: {protein_properties['mass_kda']:.2f} kDa  | pI: {protein_properties['pi']:.2f}\n"
            )
            f.write(
                f"Ext.C (Reduced):  {protein_properties['ext_coeff_reduced']} M⁻¹·cm⁻¹\n"
            )
            f.write(
                f"Ext.C (Oxidized): {protein_properties['ext_coeff_oxidized']} M⁻¹·cm⁻¹\n\n"
            )

            f.write("Sequence:\n")
            formatted_protein_sequence: str = format_sequence_ncbi_style(
                protein_record["protein_1l"]
            )
            f.write(f"{formatted_protein_sequence}\n\n\n")


def export_orf_report(
    orf_records: list[dict[str, Any]], report_file_path: str
) -> None:
    with open(report_file_path, "w", encoding="utf-8") as f:
        f.write("=" * 60 + "\n")
        f.write("SeqProfiler: ORF Analysis Report\n")
        f.write("=" * 60 + "\n\n")

        if not orf_records:
            f.write("No ORFs found matching the criteria.\n")
            return

        for orf in orf_records:
            f.write(
                f"> {orf['sequence_id']} | Strand: {orf['strand']} | Pos: {orf['start_position']} - {orf['end_position']} bp\n"
            )
            f.write("-" * 60 + "\n")

            dna_properties: dict[str, Any] = orf["dna_props"]
            f.write(
                f"DNA     Length: {dna_properties['length']} bp | Mass (ss): {dna_properties['mass_ss_da']:,.0f} Da | Mass (ds): {dna_properties['mass_ds_da']:,.0f} Da\n"
            )
            f.write(
                f"        AT: {dna_properties['at_prop']:.1f}% | GC: {dna_properties['gc_prop']:.1f}% | Tm (50mM Na+): {dna_properties['tm']:.1f} °C\n"
            )

            rna_properties: dict[str, Any] = orf["rna_props"]
            f.write(f"RNA     Mass: {rna_properties['mass_da']:,.0f} Da\n")

            protein_properties: dict[str, Any] = orf["prot_props"]
            protein_length: int = len(orf["protein_1l"])
            f.write(
                f"Protein Length: {protein_length} aa | Mass: {protein_properties['mass_kda']:.2f} kDa | pI: {protein_properties['pi']:.2f}\n"
            )
            f.write(
                f"        Ext.C (Red/Ox): {protein_properties['ext_coeff_reduced']} / {protein_properties['ext_coeff_oxidized']} M⁻¹·cm⁻¹\n\n"
            )

            f.write("Sequence:\n")
            formatted_protein_sequence: str = format_sequence_ncbi_style(orf["protein_1l"])
            f.write(f"{formatted_protein_sequence}\n\n\n\n")
