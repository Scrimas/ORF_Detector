<h1 align="center">
  SeqProfiler
</h1>

## Overview

**SeqProfiler** is a foundational bioinformatics tool written in Python. It parses FASTA files, identifies Open Reading Frames (ORFs), translates them, and calculates key biochemical properties (ssDNA/dsDNA Mass, GC%, salt-corrected Tm, pI, and dual Extinction Coefficients). It supports direct analysis of protein sequences and intelligently adapts its analysis based on the input type (e.g., mRNA vs. Genomic DNA).

## The "Vanilla Python" Philosophy

If you are reviewing this code, you might wonder why I manually hardcoded a 64-codon dictionary, built custom string-slicing loops or wrote binary search algorithms to calculate biochemical properties instead of simply importing `Bio.SeqUtils` from Biopython.

**This is intentional:** As a Master's student in Molecular and Cellular Biology (MCB) in Grenoble, my goal was to deeply understand and transcribe into code the mechanical logic of biology and the mathematics behind the biochemical properties.

## Core Features

- **Context-Aware ORF Detection:** Intelligently scans strands based on sequence type. For example, it automatically skips reverse-strand analysis for mRNA/cDNA (`NM_` accessions), as the reverse complement has no biological existence in that context.
- **Bi-Directional Scanning:** Full forward and reverse strand analysis for genomic DNA sequences.
- **Transcription & Translation:** Simulates the biological flow of information from DNA to functional protein sequences.
- **Protein Analysis:** Direct support for analyzing protein sequences (e.g., from NCBI).
- **Biochemical Calculations:** Calculates key physical properties with high biological accuracy (ssDNA/dsDNA Mass, GC%, salt-corrected Tm, pI, and both Reduced/Oxidized Extinction Coefficients).
- **NCBI Integration:** Fetch and analyze DNA, RNA, or Protein sequences directly using NCBI accession IDs.
- **Detailed Reports:** Generates structured reports presenting results in a classic bioinformatics format.

## Calculations accuracy (SeqProfiler vs. Biopython)

SeqProfiler was created as an educational exercise to translate core biological concepts into code from scratch and was certainly not made with the intention nor the pretension to outperform or replace Biopython, which remains the industry gold-standard.
And because any custom-built bioinformatics tool requires rigorous verifications, SeqProfiler includes a `pytest` suite to validate its own "Vanilla Python" algorithms against Biopython.

### Identical Logic

The following calculations are mathematically identical to Biopython's outputs:

- **GC Percentage:** Standard calculation ignoring unknown bases.
- **Extinction Coefficient:** SeqProfiler calculates both the **Reduced** and **Oxidized** extinction coefficients. The oxidized version counts 125 M⁻¹·cm⁻¹ per **pair** of cysteines (disulfide bonds), aligning exactly with Biopython's standard for proteins in oxidative environments.
- **Initiator Methionine:** Regardless of the DNA start codon used (e.g., GTG, CTG), SeqProfiler correctly translates the first amino acid of an ORF as **Methionine (M)**, reflecting the biological reality of the initiator tRNA.
- **Melting Temperature (Tm):** Matches the salt-corrected formula (assuming **50mM Na+**), providing physiologically relevant outputs for small/medium sequences.

### Intentional Divergences

While the outputs are highly accurate, there are a few intentional, well-documented biochemical divergences based on different structural assumptions:

- **Sequence Mass:** SeqProfiler reports both **Single-Stranded (ssDNA)** and **Double-Stranded (dsDNA)** mass. Biopython's default behavior assumes a 5'-phosphate group, while SeqProfiler effectively assumes a 5'-hydroxyl group. This results in a consistent 79.98 Da difference (per strand) when compared to Biopython.
- **Isoelectric Point:** Computational pI depends heavily on the dataset of pKa used. Biopython defaults to the Bjellqvist dataset whereas SeqProfiler utilizes the Lehninger pKa tables. This results in slight variations that increase with sequence length.
- **Unknown Characters:** SeqProfiler handles unknown amino acids (`X`) and nucleotides (`N`) by assigning them average biological masses (~110 Da and ~309 Da respectively) instead of ignoring them.
- **Reverse Strand Coordinates:** To respect biological directionality (5' -> 3'), ORFs found on the reverse strand are reported with a `start_position` higher than the `end_position`, mapped to the forward reference strand.

## Quick Start (Installation & Usage)

```bash
git clone https://github.com/Scrimas/SeqProfiler
cd SeqProfiler
pip install requests
```

```bash
# Using default settings
python src/main.py
```

### Available Options

| Argument         | Description                                                          | Default     |
| :--------------- | :------------------------------------------------------------------- | :---------- |
| `--min-length`   | Minimum ORF size in amino acids                                      | `50`        |
| `--input`        | Path to directory containing `.fasta` files                          | `data/`     |
| `--output`       | Path to directory for analysis reports                               | `results/`  |
| `--workers`      | Number of parallel processes to use                                  | `CPU count` |
| `--start-codons` | Comma-separated list of alternative start codons (e.g., ATG,CTG,GTG) | `ATG`       |
| `--ncbi`         | Comma-separated list of NCBI accession IDs to fetch and analyze      | `None`      |

### Examples

```bash
# Analyze local files
python src/main.py --input ./my_data
```

```bash
# Analyze sequences from NCBI
python src/main.py --ncbi NM_001301717,NC_000913
```

```bash
# Mix local and NCBI sequences
python src/main.py --input ./data --ncbi NM_001301717
```

### Running the Tests

While the test results are already included on this repository (see [tests/pytest_results.txt](https://github.com/Scrimas/SeqProfiler/blob/main/tests/pytest_results.txt)), you can verify the algorithms on your own machine.

For this, simply install `pytest` and `biopython`, then run the test suite from the root directory:

```bash
pip install pytest biopython
pytest -v
```

## Architecture

The directory structure intentionally mirrors the biological flow of information:

```text
SeqProfiler/
├── data/                       # Input .fasta files
├── results/                    # Output .txt reports
├── src/
│   ├── dna_to_codon.py         # Frame reading & Start/Stop identification
│   ├── dna_to_protein.py       # Translation logic
│   ├── dna_to_rna.py           # Transcription logic
│   ├── fasta_to_dna.py         # Sequence extraction
│   ├── main.py
│   ├── ncbi_fetch.py           # NCBI E-utilities fetching logic
│   ├── results_export.py       # Report generation
│   └── sequence_properties.py  # Biochemical properties calculations
├── tests/                      # Pytest folder
├── LICENSE
└── README.md
```

## Requirements & Compatibility

- **Python:** Python 3.14+
- **Operating Systems:**
  - **Windows 11:** Tested
  - **Linux:** Tested
  - **macOS:** Untested but should work fine as the script uses standard cross-platform libraries

## Project Status

**Completed:** This repository was created as a foundational learning exercise and portfolio piece. It is fully functional as a command-line tool and is not under strict active maintenance. It may or may not receive updates in the near or far future.

## License

This project is licensed under the MIT License.
