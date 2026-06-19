<h1 align="center">SeqProfiler</h1>

A Python tool for analyzing DNA, RNA, and protein sequences — built from scratch as part of my MCB master's degree.

SeqProfiler parses FASTA files, identifies ORFs, translates them, and computes key biochemical properties: GC%, salt-corrected Tm, ssDNA/dsDNA mass, pI, and extinction coefficients. It also fetches sequences directly from NCBI.

## Why "Vanilla Python"?

Rather than wrapping Biopython, I wrote everything by hand — codon tables, string-slicing ORF scanners, binary search pI calculations, the whole thing. The goal was to actually understand the biology and math, not just call library functions. Think of it as a translation of biology fundamentals into code.

## Features

- **Context-aware ORF detection** — skips reverse-strand analysis for mRNA/cDNA (e.g., `NM_` accessions), since the reverse complement has no biological existence in that context
- **Bi-directional scanning** for genomic DNA (forward + reverse strand)
- **Transcription & translation** following the canonical DNA → RNA → protein flow
- **Direct protein input** for sequences fetched from NCBI
- **Biochemical properties**: GC%, salt-corrected Tm (50 mM Na⁺), ssDNA/dsDNA mass, pI (Lehninger pKa tables), reduced/oxidized extinction coefficients
- **NCBI integration** via E-utilities — fetch and analyze by accession ID
- **Structured reports** in classic bioinformatics format

## Accuracy vs. Biopython

SeqProfiler ships with a `pytest` suite that benchmarks its algorithms against Biopython. Most outputs are mathematically identical; intentional divergences are documented (different pKa dataset for pI, 5'-OH instead of 5'-phosphate for mass, etc.). See [`tests/NOTES.md`](https://github.com/Scrimas/SeqProfiler/blob/main/tests/NOTES.md) for the full comparison.

## Quick Start

```bash
git clone https://github.com/Scrimas/SeqProfiler
cd SeqProfiler
pip install requests
python src/main.py
```

### Options

| Argument         | Description                                    | Default    |
| :--------------- | :--------------------------------------------- | :--------- |
| `--min-length`   | Minimum ORF size (amino acids)                 | `50`       |
| `--input`        | Directory containing `.fasta` files            | `data/`    |
| `--output`       | Directory for output reports                   | `results/` |
| `--workers`      | Parallel processes                             | CPU count  |
| `--start-codons` | Alternative start codons (e.g., `ATG,CTG,GTG`) | `ATG`      |
| `--ncbi`         | NCBI accession IDs to fetch                    | —          |

```bash
# Analyze local files
python src/main.py --input ./my_data

# Fetch from NCBI
python src/main.py --ncbi NM_001301717,NC_000913

# Mix both
python src/main.py --input ./data --ncbi NM_001301717
```

### Tests

```bash
pip install pytest biopython
pytest -v
```

## Architecture

The directory structure mirrors the biological flow of information:

```text
SeqProfiler/
├── data/                       # Input .fasta files
├── results/                    # Output reports
├── src/
│   ├── dna_to_codon.py         # Frame reading & start/stop identification
│   ├── dna_to_protein.py       # Translation
│   ├── dna_to_rna.py           # Transcription
│   ├── fasta_to_dna.py         # Sequence extraction
│   ├── main.py
│   ├── ncbi_fetch.py           # NCBI E-utilities
│   ├── results_export.py       # Report generation
│   └── sequence_properties.py  # Biochemical calculations
└── tests/
```

## Requirements

- Python 3.14+
- `requests` (runtime) — `pytest` + `biopython` for tests only
- Tested on Windows 11 and Linux

## Status

Finished and functional. Built as a learning project and portfolio piece — not under active maintenance, but may see occasional updates.

## License

MIT
