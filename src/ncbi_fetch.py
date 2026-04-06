from __future__ import annotations
import requests
import sys

def fetch_ncbi_sequence(accession_id: str) -> str:
    """
    Fetches a sequence from NCBI in FASTA format using E-utilities.
    Automatically detects if the ID belongs to 'nuccore' or 'protein' database.
    
    Args:
        accession_id (str): NCBI accession number (e.g., 'NM_001301717', 'NP_001288646').

    Returns:
        str: FASTA formatted sequence as a string.
    """
    # Simple heuristic for database selection
    # Nucleotide: NM_, NR_, NC_, NG_, XM_, XR_
    # Protein: NP_, XP_, YP_, WP_ or single letters
    
    db = "nuccore"
    protein_prefixes = ("NP_", "XP_", "YP_", "WP_", "AP_", "IP_", "ZP_")
    if accession_id.upper().startswith(protein_prefixes):
        db = "protein"
    
    url: str = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params: dict[str, str] = {
        "db": db,
        "id": accession_id,
        "rettype": "fasta",
        "retmode": "text"
    }
    
    response: requests.Response = requests.get(url, params=params)
    response.raise_for_status()
    
    fasta_data: str = response.text.strip()
    
    if not fasta_data.startswith(">"):
        raise ValueError(f"Fetched data for ID '{accession_id}' is not in FASTA format.")
        
    return fasta_data
