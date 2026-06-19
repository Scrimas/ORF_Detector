from __future__ import annotations

import requests


def fetch_ncbi_sequence(accession_number: str) -> str:
    target_database = "nuccore"
    protein_accession_prefixes = ("NP_", "XP_", "YP_", "WP_", "AP_", "IP_", "ZP_")
    if accession_number.upper().startswith(protein_accession_prefixes):
        target_database = "protein"

    eutils_url: str = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    request_parameters: dict[str, str] = {
        "db": target_database,
        "id": accession_number,
        "rettype": "fasta",
        "retmode": "text",
    }

    http_response: requests.Response = requests.get(eutils_url, params=request_parameters)
    http_response.raise_for_status()

    fetched_fasta: str = http_response.text.strip()

    if not fetched_fasta.startswith(">"):
        raise ValueError(
            f"Fetched data for ID '{accession_number}' is not in FASTA format."
        )

    return fetched_fasta
