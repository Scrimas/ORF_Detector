## Overview

**ORF_Detector** is a purely foundational, from-scratch bioinformatics tool written in Python. It parses standard FASTA files, identifies valid Open Reading Frames (ORFs), transcribes DNA into mRNA, and translates sequences into both 3-letter and 1-letter amino acid structures.

## The "Vanilla Python" Philosophy

If you are reviewing this code, you might wonder why I manually hardcoded a 64-codon dictionary or built custom string-slicing loops to map reading frames instead of simply importing `Bio.Seq` from Biopython.

**This is an intentional and pedagogical constraint.** As a 3rd-year biology student in Grenoble, my goal was to deeply understand and computationally replicate the mechanical logic of the ribosomal machinery. Relying on "black box" libraries would defeat the purpose of the exercise. Every algorithmic choice (from the step-of-3 indexing to the reverse-complement strand mathematics) was written explicitly to prove a granular understanding of genomic translation.

## Core Features

- **Custom Multi-FASTA Parsing:** Efficiently handles multi-line sequence concatenation and whitespace stripping without external dependencies.
    
- **Strict ORF Detection:** Scans both forward (+) and reverse (-) strands for the `ATG` start codon, stepping exactly in-frame (modulo 3) until a stop codon (`TAA`, `TAG`, `TGA`) is reached.
    
- **Intelligent Redundancy Filtering:** Utilizes a heuristic `used_stops` tracking system to prioritize the longest functional reading frame and prevent reporting smaller, nested ORFs.
    
- **Bi-Directional Coordinate Mapping:** Accurately maps the 1-based start/end coordinates of reverse-complement ORFs back to the original forward strand reference.
    
- **Batch Processing:** Automatically processes all `.fasta` files in the `/data` directory and exports structured reports to `/results`.

## Architecture

The directory structure intentionally mirrors the biological flow of information:

```
ORF_Detector/
├── data/                  # Input .fasta files
├── results/               # Output text reports
└── src/
    ├── fasta_to_dna.py    # Sequence extraction & multi-FASTA mapping
    ├── dna_to_codon.py    # Frame reading & Start/Stop identification
    ├── dna_to_rna.py      # Transcription logic
    ├── dna_to_protein.py  # Translation (1-letter & 3-letter formats)
    ├── results_export.py  # Structured report generation
    └── main.py            # The central execution orchestrator
```

## Usage

1. **Prepare Data:** Place your target `.fasta` files into the `data/` directory.
    
2. **Run Pipeline:** Execute the main script from the root directory:

    ```
    python src/main.py
    ```

3. **Configure:** Enter your desired minimum ORF size in amino acids (defaults to 50) when prompted.

4. **Review:** Retrieve detailed analysis reports in the `results/` folder.

## About the Author

**Ismaël PHILIPPE**

_3rd-year Biology Student | Grenoble, France_

I am focused on bridging the gap between life sciences and computer science. This project serves as a practical demonstration of translating complex biological paradigms into functional, modular code.


## License

This project is licensed under the MIT License.