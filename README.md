<h1 align="center">
  SeqProfiler
</h1>

## Overview

**SeqProfiler** is a purely foundational, from-scratch bioinformatics tool written in Python. It parses standard FASTA files, identifies valid Open Reading Frames (ORFs), transcribes DNA into mRNA, translates sequences into amino acids, and calculates key biochemical properties (Mass, Tm, pI, Extinction Coefficient) for all identified molecules.

## Project Status

**Completed / Open to updates:** This repository was created as a foundational learning exercise and portfolio piece. While it is fully functional as a command-line tool and is not under strict active maintenance, a potential future update I am exploring is wrapping the pipeline in a Graphical User Interface (GUI).

## The "Vanilla Python" Philosophy

If you are reviewing this code, you might wonder why I manually hardcoded a 64-codon dictionary, built custom string-slicing loops, or wrote binary search algorithms to calculate isoelectric points instead of simply importing `Bio.Seq` from Biopython.

**This is an intentional and pedagogical constraint.** As a 3rd-year biology student in Grenoble, my goal was to deeply understand and computationally replicate the mechanical logic of the ribosomal machinery and the mathematics behind biochemical properties. Relying on "black box" libraries would defeat the purpose of the exercise. Every algorithmic choice was written explicitly to prove a granular understanding of genomic translation and molecular biophysics.

## Core Features

- **Batch Sequence Analysis:** Simultaneously processes multiple FASTA files using parallel computing to handle large-scale genomic data efficiently.
- **Bi-Directional ORF Detection:** Automatically scans both forward and reverse DNA strands to identify all potential protein-coding regions (ORFs).
- **Automated Transcription & Translation:** Converts identified DNA sequences into mRNA and provides protein translations in both 1-letter and 3-letter formats.
- **Biochemical Characterization:** Calculates key physical properties including molecular mass, melting temperature (Tm), isoelectric point (pI), and extinction coefficient.
- **Detailed Reporting:** Generates structured analysis reports for every identified molecule, mapping ORF coordinates back to the original reference sequence and presenting protein sequences in NCBI-style chunked blocks with position numbers.

## Architecture

The directory structure intentionally mirrors the biological flow of information, now enhanced with modern Python type hinting and optimized algorithms:

```text
SeqProfiler/
├── data/                       # Input .fasta files
├── results/                    # Output text reports
└── src/
    ├── dna_to_codon.py         # Frame reading & Start/Stop identification
    ├── dna_to_protein.py       # Translation (1-letter & 3-letter formats)
    ├── dna_to_rna.py           # Transcription logic
    ├── fasta_to_dna.py         # Sequence extraction & mapping
    ├── main.py                 # Main execution point
    ├── results_export.py       # Structured & NCBI-style report generation
    └── sequence_properties.py  # Biochemical property calculations
```

## Usage

1. **Clone repository**

```bash
git clone https://github.com/Scrimas/SeqProfiler
```
    
2. **Run Pipeline:** Execute the main script from the root directory:

```bash
# Using default settings
python src/main.py
```

```bash
# Specifying custom settings
python src/main.py --min-length 100 --workers 4 --input ./data_folder --output ./results_folder 
```

3. **Review:** Retrieve detailed analysis reports in your specified output folder (defaults to `results/`).

### Available Options

| Argument | Description | Default |
| :--- | :--- | :--- |
| `--min-length` | Minimum ORF size in amino acids | `50` |
| `--input` | Path to directory containing `.fasta` files | `data/` |
| `--output` | Path to directory for analysis reports | `results/` |
| `--workers` | Number of parallel processes to use | `CPU count` |
| `--start-codons` | Comma-separated list of alternative start codons (e.g., ATG,CTG,GTG) | `ATG` |

## License

This project is licensed under the MIT License.
