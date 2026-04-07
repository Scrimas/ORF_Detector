import pytest
from sequence_properties import (
    calculate_dna_properties,
    calculate_rna_properties,
    calculate_protein_properties,
)
from dna_to_protein import translate_sequence
from Bio.SeqUtils import molecular_weight, gc_fraction
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# --- TEST DATA ---
SHORT_DNA = "ATGTGGTACTGCTGCGACGAAAAACGCCACGGCGCC"
SHORT_RNA = SHORT_DNA.replace("T", "U")
SHORT_PROT = translate_sequence(SHORT_DNA)[1].replace("*", "")

LONG_DNA = "GTTCAGTGGTGGCTTTTGGCCAAGTGCCCGCAACCGAGCCTCGGCTCCGCCCGGCAGTTAGTAAGCTCCTTCGGAGGGGACAGTGGTAAAGCCGACTTGATCGGGACAGACTTCGAGAATGGTTAACCGTAGAGTAGTTACATAAAGGATATTCAGAGGTGGTGATGTCACCTAAAGGCCTACAAGATTAGCCTGCTTTGCTCTCAGATGACCATCTGTGATCACTTACGATAATGCAGGGCGCAGTACGATGACCTCATATGACGATGGGTCGGTACTAATAAACGAATAGAAGAGTGCCGCCTCGCCTATCTGTATCGGAATTAGCTCGTGAAAGCCGTATCACAACCGACCCTGGACGATATTTGGCTAGTCGGTCTTCTTAACGTCATCGCTAATTACCCCACTCCTCGAGGACGTCCAGTCACAAGCCACTCTTTAGCGTCTCGGGTCAGATCGTTGGTTGATCTTGGGGAGCTATTATGGATCACGCTATTTAAGGCAGGAGCGTAATACGCATTTCTTGGCTATAGTACCACAGGTAGCCGACGTAGGCACAAGGCCGAGTAACGGGAGAGAATTTGCGTCTCCCCATGTCTAATTATATAGTTTAAGAATAACGCATACCCCCCTGGGTGCAACGGGCTTGTGATGTCTTACGTGTTGTGGGCACGGGTGAGAACCCTAGTACCCATTGTATATTCCGAATCAAACCGGGAGCCTCTATCTCGTACCATCCTGCGATCACCGCCTCACTAGCCCACACCTTCCCGCGTGGCCAAATACAGACACATGGCAGGTAGAGTAACAAAGAGAGGTATCGCGGTGCAGGTTTTACTGTTGTACTGTACAACCGGCTAAGCGAAAGCGCGGGTTTGGCTTCGGTGTAGCCCTCGGTAGAGCCTGGCATAGGCTTTGAGCGGGCGGCTATCCTTAGGGCTCGCACTGTATTTCATTAGTACGGACTTGTCAGGGATGTAATCACGTCCGAGTAGTGAGCAG"
LONG_RNA = LONG_DNA.replace("T", "U")
LONG_PROT = translate_sequence(LONG_DNA)[1].replace("*", "")


class TestDNAProperties:
    """Groups all DNA-related tests together."""

    def test_short_dna_mass(self):
        user_props = calculate_dna_properties(SHORT_DNA)
        bio_mass = molecular_weight(
            SHORT_DNA, "DNA", circular=False, double_stranded=False
        )
        assert user_props["mass_da"] == pytest.approx(bio_mass - 79.98, rel=1e-3)

    def test_short_dna_gc_content(self):
        user_props = calculate_dna_properties(SHORT_DNA)
        bio_gc = gc_fraction(SHORT_DNA) * 100
        assert user_props["gc_prop"] == pytest.approx(bio_gc, rel=1e-4)

    def test_short_dna_melting_temp(self):
        user_props = calculate_dna_properties(SHORT_DNA)
        bio_tm = 64.9 + 41 * (SHORT_DNA.count("G") + SHORT_DNA.count("C") - 16.4) / len(
            SHORT_DNA
        )
        assert user_props["tm"] == pytest.approx(bio_tm, rel=1e-4)

    def test_long_dna_mass(self):
        user_props = calculate_dna_properties(LONG_DNA)
        bio_mass = molecular_weight(
            LONG_DNA, "DNA", circular=False, double_stranded=False
        )
        assert user_props["mass_da"] == pytest.approx(bio_mass - 79.98, rel=1e-3)

    def test_long_dna_gc_content(self):
        user_props = calculate_dna_properties(LONG_DNA)
        bio_gc = gc_fraction(LONG_DNA) * 100
        assert user_props["gc_prop"] == pytest.approx(bio_gc, rel=1e-4)

    def test_long_dna_melting_temp(self):
        user_props = calculate_dna_properties(LONG_DNA)
        bio_tm = 64.9 + 41 * (LONG_DNA.count("G") + LONG_DNA.count("C") - 16.4) / len(
            LONG_DNA
        )
        assert user_props["tm"] == pytest.approx(bio_tm, rel=1e-4)


class TestRNAProperties:
    """Groups all RNA-related tests together."""

    def test_short_rna_mass(self):
        user_props = calculate_rna_properties(SHORT_RNA)
        bio_mass = molecular_weight(
            SHORT_RNA, "RNA", circular=False, double_stranded=False
        )
        assert user_props["mass_da"] == pytest.approx(bio_mass - 79.98, rel=1e-3)

    def test_long_rna_mass(self):
        user_props = calculate_rna_properties(LONG_RNA)
        bio_mass = molecular_weight(
            LONG_RNA, "RNA", circular=False, double_stranded=False
        )
        assert user_props["mass_da"] == pytest.approx(bio_mass - 79.98, rel=1e-3)


class TestProteinProperties:
    """Groups all Protein-related tests together."""

    def test_short_protein_mass(self):
        user_props = calculate_protein_properties(SHORT_PROT)
        bio_mass = molecular_weight(SHORT_PROT, "protein", circular=False) / 1000.0
        assert user_props["mass_kda"] == pytest.approx(bio_mass, rel=1e-3)

    def test_short_protein_pi(self):
        user_props = calculate_protein_properties(SHORT_PROT)
        bio_pi = ProteinAnalysis(SHORT_PROT).isoelectric_point()
        assert user_props["pi"] == pytest.approx(bio_pi, abs=0.9)

    def test_short_protein_extinction_coefficient(self):
        user_props = calculate_protein_properties(SHORT_PROT)
        _, bio_ext_oxidized = ProteinAnalysis(SHORT_PROT).molar_extinction_coefficient()

        cysteine_count = SHORT_PROT.count("C")
        biopython_cysteine_addition = (cysteine_count // 2) * 125
        seqprofiler_cysteine_addition = cysteine_count * 125
        expected_difference = (
            seqprofiler_cysteine_addition - biopython_cysteine_addition
        )

        assert user_props["ext_coeff"] == pytest.approx(
            bio_ext_oxidized + expected_difference, rel=1e-3
        )

    def test_long_protein_mass(self):
        user_props = calculate_protein_properties(LONG_PROT)
        bio_mass = molecular_weight(LONG_PROT, "protein", circular=False) / 1000.0
        assert user_props["mass_kda"] == pytest.approx(bio_mass, rel=1e-3)

    def test_long_protein_pi(self):
        user_props = calculate_protein_properties(LONG_PROT)
        bio_pi = ProteinAnalysis(LONG_PROT).isoelectric_point()

        assert user_props["pi"] == pytest.approx(bio_pi, abs=0.4)

    def test_long_protein_extinction_coefficient(self):
        user_props = calculate_protein_properties(LONG_PROT)
        _, bio_ext_oxidized = ProteinAnalysis(LONG_PROT).molar_extinction_coefficient()

        cysteine_count = LONG_PROT.count("C")
        biopython_cysteine_addition = (cysteine_count // 2) * 125
        seqprofiler_cysteine_addition = cysteine_count * 125
        expected_difference = (
            seqprofiler_cysteine_addition - biopython_cysteine_addition
        )

        assert user_props["ext_coeff"] == pytest.approx(
            bio_ext_oxidized + expected_difference, rel=1e-3
        )
