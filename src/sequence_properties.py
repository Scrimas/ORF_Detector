from __future__ import annotations


def compute_dna_properties(dna_sequence: str) -> dict[str, int | float]:
    sequence_length: int = len(dna_sequence)
    if sequence_length == 0:
        return {
            "length": 0,
            "mass_ss_da": 0.0,
            "mass_ds_da": 0.0,
            "gc_prop": 0.0,
            "at_prop": 0.0,
            "tm": 0.0,
        }

    count_a: int = dna_sequence.count("A")
    count_t: int = dna_sequence.count("T")
    count_c: int = dna_sequence.count("C")
    count_g: int = dna_sequence.count("G")
    count_n: int = dna_sequence.count("N")

    known_base_count: int = count_a + count_t + count_c + count_g
    if known_base_count > 0:
        gc_percentage: float = ((count_g + count_c) / known_base_count) * 100
        at_percentage: float = ((count_a + count_t) / known_base_count) * 100
    else:
        gc_percentage, at_percentage = 0.0, 0.0

    single_stranded_mass: float = (
        (count_a * 313.21)
        + (count_t * 304.2)
        + (count_c * 289.18)
        + (count_g * 329.21)
        + (count_n * 308.95)
        - 61.96
    )
    double_stranded_mass: float = (
        single_stranded_mass
        + (count_t * 313.21)
        + (count_a * 304.2)
        + (count_g * 289.18)
        + (count_c * 329.21)
        + (count_n * 308.95)
        - 61.96
    )

    import math

    sodium_concentration = 0.05
    melting_temperature: float = (
        81.5
        + 16.6 * math.log10(sodium_concentration)
        + 0.41 * gc_percentage
        - 675 / sequence_length
    )

    return {
        "length": sequence_length,
        "mass_ss_da": single_stranded_mass,
        "mass_ds_da": double_stranded_mass,
        "gc_prop": gc_percentage,
        "at_prop": at_percentage,
        "tm": melting_temperature,
    }


def compute_rna_properties(rna_sequence: str) -> dict[str, int | float]:
    sequence_length: int = len(rna_sequence)
    if sequence_length == 0:
        return {"length": 0, "mass_da": 0.0}

    count_a: int = rna_sequence.count("A")
    count_u: int = rna_sequence.count("U")
    count_c: int = rna_sequence.count("C")
    count_g: int = rna_sequence.count("G")
    count_n: int = rna_sequence.count("N")

    molecular_weight: float = (
        (count_a * 329.21)
        + (count_u * 306.15)
        + (count_c * 305.18)
        + (count_g * 345.21)
        + (count_n * 321.44)
        - 61.96
    )
    return {"length": sequence_length, "mass_da": molecular_weight}


def compute_protein_properties(protein_sequence: str) -> dict[str, float]:
    clean_protein_sequence: str = protein_sequence.replace("*", "").replace("-", "")
    sequence_length: int = len(clean_protein_sequence)
    if sequence_length == 0:
        return {
            "mass_kda": 0.0,
            "pi": 0.0,
            "ext_coeff_reduced": 0.0,
            "ext_coeff_oxidized": 0.0,
        }

    amino_acid_mass_table: dict[str, float] = {
        "A": 71.0788,
        "R": 156.1875,
        "N": 114.1038,
        "D": 115.0886,
        "C": 103.1388,
        "E": 129.1155,
        "Q": 128.1307,
        "G": 57.0519,
        "H": 137.1411,
        "I": 113.1594,
        "L": 113.1594,
        "K": 128.1741,
        "M": 131.1926,
        "F": 147.1766,
        "P": 97.1167,
        "S": 87.0782,
        "T": 101.1051,
        "W": 186.2132,
        "Y": 163.1760,
        "V": 99.1326,
        "X": 110.0,
    }

    molecular_weight_da: float = (
        sum(amino_acid_mass_table.get(aa, 110.0) for aa in clean_protein_sequence)
        + 18.01524
    )
    molecular_weight_kda: float = molecular_weight_da / 1000.0

    amino_acid_counts: dict[str, int] = {
        aa: clean_protein_sequence.count(aa) for aa in set(clean_protein_sequence)
    }

    tryptophan_count: int = amino_acid_counts.get("W", 0)
    tyrosine_count: int = amino_acid_counts.get("Y", 0)
    cysteine_count: int = amino_acid_counts.get("C", 0)

    reduced_extinction_coefficient: float = (tryptophan_count * 5500) + (
        tyrosine_count * 1490
    )
    oxidized_extinction_coefficient: float = reduced_extinction_coefficient + (
        (cysteine_count // 2) * 125
    )

    pka_n_terminus: float = 9.69
    pka_c_terminus: float = 2.34
    basic_pka_table: dict[str, float] = {"K": 10.53, "R": 12.48, "H": 6.00}
    acidic_pka_table: dict[str, float] = {"D": 3.65, "E": 4.25, "C": 8.18, "Y": 10.07}

    def calculate_net_charge(ph_value: float) -> float:
        calculated_charge: float = 0.0

        calculated_charge += (
            10 ** (pka_n_terminus - ph_value) / (1 + 10 ** (pka_n_terminus - ph_value))
        )
        for aa, pka in basic_pka_table.items():
            calculated_charge += amino_acid_counts.get(aa, 0) * (
                10 ** (pka - ph_value) / (1 + 10 ** (pka - ph_value))
            )

        calculated_charge -= (
            10 ** (ph_value - pka_c_terminus) / (1 + 10 ** (ph_value - pka_c_terminus))
        )
        for aa, pka in acidic_pka_table.items():
            calculated_charge -= amino_acid_counts.get(aa, 0) * (
                10 ** (ph_value - pka) / (1 + 10 ** (ph_value - pka))
            )
        return calculated_charge

    lower_bound: float = 0.0
    upper_bound: float = 14.0
    isoelectric_point: float = 7.0
    for _ in range(100):
        isoelectric_point = (lower_bound + upper_bound) / 2
        if calculate_net_charge(isoelectric_point) > 0:
            lower_bound = isoelectric_point
        else:
            upper_bound = isoelectric_point

    return {
        "mass_kda": molecular_weight_kda,
        "pi": isoelectric_point,
        "ext_coeff_reduced": reduced_extinction_coefficient,
        "ext_coeff_oxidized": oxidized_extinction_coefficient,
    }
