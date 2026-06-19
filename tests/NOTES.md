# Implementation Notes

Technical decisions and divergences worth documenting, for anyone reading the source or comparing outputs against other tools.

## Accuracy vs. Biopython

SeqProfiler's algorithms are benchmarked against Biopython in `tests/`. Here's what matches, what diverges, and why.

### Identical outputs

**GC%** — Standard calculation, unknown bases ignored. Matches exactly.

**Extinction coefficient** — Both reduced and oxidized variants are computed. The oxidized version counts 125 M⁻¹·cm⁻¹ per _pair_ of cysteines (disulfide bonds), which aligns with Biopython's convention for proteins in oxidative environments.

**Initiator methionine** — Regardless of which start codon is used (ATG, GTG, CTG…), the first amino acid of every ORF is translated as methionine. This reflects the biology: the initiator tRNA carries Met regardless of the codon it reads.

**Melting temperature (Tm)** — Uses the salt-corrected formula assuming 50 mM Na⁺. Matches Biopython for small to medium sequences.

### Intentional divergences

**Sequence mass (ssDNA / dsDNA)** —
Biopython assumes a 5'-phosphate group by default. SeqProfiler assumes a 5'-hydroxyl. This produces a consistent ~79.98 Da difference per strand. Neither is wrong — they're just different structural assumptions about the terminus.

**Isoelectric point (pI)** —
Biopython uses the Bjellqvist pKa dataset. SeqProfiler uses the Lehninger tables. Both are valid reference datasets; the gap grows with sequence length as small per-residue differences accumulate.

**Unknown characters** —
SeqProfiler assigns average biological masses to unknown residues rather than skipping them: `X` (unknown amino acid) → ~110 Da, `N` (unknown nucleotide) → ~309 Da. Biopython ignores unknowns, which can silently undercount mass for sequences with ambiguity codes.

**Reverse strand coordinates** —
ORFs on the reverse strand are reported with `start_position > end_position`, mapped to the forward reference strand. This preserves biological directionality (5' → 3') rather than always returning a lower start than end, which would lose orientation information.
