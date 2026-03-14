def calculate_dna_properties(sequence):
    length = len(sequence)
    if length == 0:
        return {"length": 0, "mass_da": 0, "gc_prop": 0, "at_prop": 0, "tm": 0}
    a = sequence.count('A')
    t = sequence.count('T')
    c = sequence.count('C')
    g = sequence.count('G')
    gc_prop = ((g + c) / length) * 100
    at_prop = ((a + t) / length) * 100
    mass = (a * 313.21) + (t * 304.2) + (c * 289.18) + (g * 329.21) - 61.96
    if length < 14:
        tm = (a + t) * 2 + (g + c) * 4
    else:
        tm = 64.9 + 41 * (g + c - 16.4) / length
    return {"length": length, "mass_da": mass, "gc_prop": gc_prop, "at_prop": at_prop, "tm": tm}

def calculate_rna_properties(sequence):
    length = len(sequence)
    if length == 0:
        return {"length": 0, "mass_da": 0}
    a = sequence.count('A')
    u = sequence.count('U')
    c = sequence.count('C')
    g = sequence.count('G')
    mass = (a * 329.21) + (u * 306.15) + (c * 305.18) + (g * 345.21) - 61.96
    return {"length": length, "mass_da": mass}

def calculate_protein_properties(sequence):
    seq = sequence.replace('*', '')
    length = len(seq)
    if length == 0:
        return {"mass_kda": 0, "pi": 0, "ext_coeff": 0}
    aa_masses = {
        'A': 71.0788, 'R': 156.1875, 'N': 114.1038, 'D': 115.0886,
        'C': 103.1388, 'E': 129.1155, 'Q': 128.1307, 'G': 57.0519,
        'H': 137.1411, 'I': 113.1594, 'L': 113.1594, 'K': 128.1741,
        'M': 131.1926, 'F': 147.1766, 'P': 97.1167, 'S': 87.0782,
        'T': 101.1051, 'W': 186.2132, 'Y': 163.1760, 'V': 99.1326, 'X': 0
    }
    mass_da = sum(aa_masses.get(aa, 0) for aa in seq) + 18.01524
    mass_kda = mass_da / 1000.0
    w = seq.count('W')
    y = seq.count('Y')
    c = seq.count('C')
    ext_coeff = (w * 5500) + (y * 1490) + (c * 125)
    pka_n_term = 9.69
    pka_c_term = 2.34
    pka_basic = {'K': 10.53, 'R': 12.48, 'H': 6.00}
    pka_acidic = {'D': 3.65, 'E': 4.25, 'C': 8.18, 'Y': 10.07}

    def net_charge(pH):
        charge = 0.0
        charge += 10**(pka_n_term - pH) / (1 + 10**(pka_n_term - pH))
        for aa, pka in pka_basic.items():
            charge += seq.count(aa) * (10**(pka - pH) / (1 + 10**(pka - pH)))
        charge -= 10**(pH - pka_c_term) / (1 + 10**(pH - pka_c_term))
        for aa, pka in pka_acidic.items():
            charge -= seq.count(aa) * (10**(pH - pka) / (1 + 10**(pH - pka)))
        return charge

    low, high = 0.0, 14.0
    pi = 7.0
    for _ in range(100):
        pi = (low + high) / 2
        if net_charge(pi) > 0:
            low = pi
        else:
            high = pi
    return {"mass_kda": mass_kda, "pi": pi, "ext_coeff": ext_coeff}