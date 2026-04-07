<h1 align="center">
  SeqProfiler
</h1>

## Overview

**SeqProfiler** is a purely foundational, from-scratch bioinformatics tool written in Python. It parses FASTA files, identifies Open Reading Frames (ORFs) in DNA/RNA, translates them into amino acids sequences, and calculates key biochemical properties (Mass, GC percentage, Tm, pI and Extinction Coefficient). It also supports direct analysis of protein sequences.

## The "Vanilla Python" Philosophy

If you are reviewing this code, you might wonder why I manually hardcoded a 64-codon dictionary, built custom string-slicing loops or wrote binary search algorithms to calculate biochemical properties instead of simply importing `Bio.SeqUtils` from Biopython.

**This is intentional:** As a 3rd-year B.Sc. biology student in Grenoble, my goal was to deeply understand, but more importantly, to transcribe into code the mechanical logic of biology and the mathematics behind the biochemical properties. Relying on "black box" libraries would defeat the entire purpose of the exercise.

## Core Features

- **Bi-Directional ORFs Detection:** Scans both forward and reverse DNA/RNA strands to identify all potential protein-coding regions.
- **Transcription & Translation:** Transcribes and translates identified ORFs into 1-letter protein sequences.
- **Protein Analysis:** Direct support for analyzing protein sequences (e.g., from NCBI).
- **Biochemical calculations:** Calculates key physical properties (Mass, GC percentage, Tm, pI and Extinction Coefficient).
- **NCBI Integration:** Fetch and analyze DNA, RNA, or Protein sequences directly using NCBI accession IDs.
- **Detailed Reports:** Generates structured reports for every ORF or Protein, maps coordinates and presents 1-letter protein sequences in a NCBI-style chunked blocks with position numbers.

## Calculations accuracy (SeqProfiler vs. Biopython)

SeqProfiler was created as an educational exercise to translate core biological concepts into code from scratch and was certainly not made with the intention nor the pretension to outperform or replace Biopython, which remains the industry gold-standard.
And because any custom-built bioinformatics tool requires rigorous verifications, SeqProfiler includes a `pytest` suite to validate its own "Vanilla Python" algorithms against Biopython to prove the fundamental math is sound.

While the outputs are highly accurate, there are a few intentional, well-documented biochemical divergences based on different structural assumptions :

- **Sequence Mass:** Biopython's default behavior assumes a 5'-phosphate group while SeqProfiler's algorithm effectively assumes a 5'-hydroxyl group.
  This results in a consistent 79.98 Da difference.

- **Isoelectric Point:** Computational pI depends heavily on the dataset of pKa used. Biopython defaults to the Bjellqvist dataset whereas SeqProfiler utilizes the Lehninger pKa tables.
  This results in slight increasing pI variations the longer the sequence is.

- **Extinction Coefficient:** Biopython defaults to calculating an "oxidized" state, meaning it only counts cysteines that form disulfide bonds, adding a coefficient of 125 M⁻¹·cm⁻¹ per pair of cysteines. SeqProfiler takes a generalized approach, adding a coefficient of 125 M⁻¹·cm⁻¹ for every individual cysteine present in the sequence.
  This results in increasing variations in extinction coefficient the more cysteines the sequence contains.

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

# Analyze sequences from NCBI
python src/main.py --ncbi NM_001301717,NC_000913

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
