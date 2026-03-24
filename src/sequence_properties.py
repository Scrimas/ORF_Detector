from __future__ import annotations

def calculate_dna_properties(sequence: str) -> dict[str, int | float]:
    """
    Calculates physical properties for a DNA sequence.
    
    Args:
        sequence (str): The DNA sequence.
        
    Returns:
        Dict: Length, mass (Da), GC%, AT%, and Melting Temp (Tm).
    """
    length: int = len(sequence)
    if length == 0:
        return {"length": 0, "mass_da": 0.0, "gc_prop": 0.0, "at_prop": 0.0, "tm": 0.0}
        
    a: int = sequence.count('A')
    t: int = sequence.count('T')
    c: int = sequence.count('C')
    g: int = sequence.count('G')
    
    gc_prop: float = ((g + c) / length) * 100
    at_prop: float = ((a + t) / length) * 100
    
    mass: float = (a * 313.21) + (t * 304.2) + (c * 289.18) + (g * 329.21) - 61.96
    
    tm: float
    if length < 14:
        tm = (a + t) * 2 + (g + c) * 4
    else:
        tm = 64.9 + 41 * (g + c - 16.4) / length
        
    return {"length": length, "mass_da": mass, "gc_prop": gc_prop, "at_prop": at_prop, "tm": tm}

def calculate_rna_properties(sequence: str) -> dict[str, int | float]:
    """
    Calculates physical properties for an mRNA sequence.
    
    Args:
        sequence (str): The RNA sequence.
        
    Returns:
        Dict: Length and mass (Da).
    """
    length: int = len(sequence)
    if length == 0:
        return {"length": 0, "mass_da": 0.0}
        
    a: int = sequence.count('A')
    u: int = sequence.count('U')
    c: int = sequence.count('C')
    g: int = sequence.count('G')
    
    mass: float = (a * 329.21) + (u * 306.15) + (c * 305.18) + (g * 345.21) - 61.96
    return {"length": length, "mass_da": mass}

def calculate_protein_properties(sequence: str) -> dict[str, float]:
    """
    Calculates biochemical properties for a protein sequence.
    
    Optimized pI calculation by pre-counting amino acids to avoid O(N) 
    operations within the binary search loop.
    
    Args:
        sequence (str): The 1-letter protein sequence.
        
    Returns:
        Dict: Mass (kDa), pI, and Extinction Coefficient.
    """

    seq: str = sequence.replace('*', '')
    length: int = len(seq)
    if length == 0:
        return {"mass_kda": 0.0, "pi": 0.0, "ext_coeff": 0.0}
        
    aa_masses: dict[str, float] = {
        'A': 71.0788, 'R': 156.1875, 'N': 114.1038, 'D': 115.0886,
        'C': 103.1388, 'E': 129.1155, 'Q': 128.1307, 'G': 57.0519,
        'H': 137.1411, 'I': 113.1594, 'L': 113.1594, 'K': 128.1741,
        'M': 131.1926, 'F': 147.1766, 'P': 97.1167, 'S': 87.0782,
        'T': 101.1051, 'W': 186.2132, 'Y': 163.1760, 'V': 99.1326, 'X': 0.0
    }
    
    # Molecular Weight calculations
    mass_da: float = sum(aa_masses.get(aa, 0.0) for aa in seq) + 18.01524
    mass_kda: float = mass_da / 1000.0
    
    aa_counts: dict[str, int] = {aa: seq.count(aa) for aa in set(seq)}
    
    # Extinction Coefficient calculations
    w: int = aa_counts.get('W', 0)
    y: int = aa_counts.get('Y', 0)
    c: int = aa_counts.get('C', 0)
    ext_coeff: float = (w * 5500) + (y * 1490) + (c * 125)
    
    # pI calculation Constants
    pka_n_term: float = 9.69
    pka_c_term: float = 2.34
    pka_basic: dict[str, float] = {'K': 10.53, 'R': 12.48, 'H': 6.00}
    pka_acidic: dict[str, float] = {'D': 3.65, 'E': 4.25, 'C': 8.18, 'Y': 10.07}

    def net_charge(pH: float) -> float:
        """Helper to calculate net charge at a given pH."""
        charge: float = 0.0

        # Positive charges (N-term + basic AAs)
        charge += 10**(pka_n_term - pH) / (1 + 10**(pka_n_term - pH))
        for aa, pka in pka_basic.items():
            charge += aa_counts.get(aa, 0) * (10**(pka - pH) / (1 + 10**(pka - pH)))

        # Negative charges (C-term + acidic AAs)
        charge -= 10**(pH - pka_c_term) / (1 + 10**(pH - pka_c_term))
        for aa, pka in pka_acidic.items():
            charge -= aa_counts.get(aa, 0) * (10**(pH - pka) / (1 + 10**(pH - pka)))
        return charge

    low: float = 0.0
    high: float = 14.0
    pi: float = 7.0
    for _ in range(100):
        pi = (low + high) / 2
        if net_charge(pi) > 0:
            low = pi
        else:
            high = pi
            
    return {"mass_kda": mass_kda, "pi": pi, "ext_coeff": ext_coeff}
