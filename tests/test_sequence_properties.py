import pytest
import math
from Bio.SeqUtils import gc_fraction, molecular_weight
from Bio.SeqUtils.ProtParam import ProteinAnalysis

from dna_to_protein import translate_sequence
from sequence_properties import (
    compute_dna_properties,
    compute_protein_properties,
    compute_rna_properties,
)

SHORT_DNA = "ATGTGGTACTGCTGCGACGAAAAACGCCACGGCGCC"
SHORT_RNA = SHORT_DNA.replace("T", "U")
SHORT_PROT = translate_sequence(SHORT_RNA)[1].replace("*", "")

LONG_DNA = "GTTCAGTGGTGGCTTTTGGCCAAGTGCCCGCAACCGAGCCTCGGCTCCGCCCGGCAGTTAGTAAGCTCCTTCGGAGGGGACAGTGGTAAAGCCGACTTGATCGGGACAGACTTCGAGAATGGTTAACCGTAGAGTAGTTACATAAAGGATATTCAGAGGTGGTGATGTCACCTAAAGGCCTACAAGATTAGCCTGCTTTGCTCTCAGATGACCATCTGTGATCACTTACGATAATGCAGGGCGCAGTACGATGACCTCATATGACGATGGGTCGGTACTAATAAACGAATAGAAGAGTGCCGCCTCGCCTATCTGTATCGGAATTAGCTCGTGAAAGCCGTATCACAACCGACCCTGGACGATATTTGGCTAGTCGGTCTTCTTAACGTCATCGCTAATTACCCCACTCCTCGAGGACGTCCAGTCACAAGCCACTCTTTAGCGTCTCGGGTCAGATCGTTGGTTGATCTTGGGGAGCTATTATGGATCACGCTATTTAAGGCAGGAGCGTAATACGCATTTCTTGGCTATAGTACCACAGGTAGCCGACGTAGGCACAAGGCCGAGTAACGGGAGAGAATTTGCGTCTCCCCATGTCTAATTATATAGTTTAAGAATAACGCATACCCCCCTGGGTGCAACGGGCTTGTGATGTCTTACGTGTTGTGGGCACGGGTGAGAACCCTAGTACCCATTGTATATTCCGAATCAAACCGGGAGCCTCTATCTCGTACCATCCTGCGATCACCGCCTCACTAGCCCACACCTTCCCGCGTGGCCAAATACAGACACATGGCAGGTAGAGTAACAAAGAGAGGTATCGCGGTGCAGGTTTTACTGTTGTACTGTACAACCGGCTAAGCGAAAGCGCGGGTTTGGCTTCGGTGTAGCCCTCGGTAGAGCCTGGCATAGGCTTTGAGCGGGCGGCTATCCTTAGGGCTCGCACTGTATTTCATTAGTACGGACTTGTCAGGGATGTAATCACGTCCGAGTAGTGAGCAG"
LONG_RNA = LONG_DNA.replace("T", "U")
LONG_PROT = translate_sequence(LONG_RNA)[1].replace("*", "")


class TestDNAProperties:

    def test_short_dna_mass(self):
        computed_properties = compute_dna_properties(SHORT_DNA)
        bio_mass_ds = molecular_weight(
            SHORT_DNA, "DNA", circular=False, double_stranded=True
        )
        bio_mass_ss = molecular_weight(
            SHORT_DNA, "DNA", circular=False, double_stranded=False
        )
        assert computed_properties["mass_ds_da"] == pytest.approx(bio_mass_ds - 159.96, rel=1e-3)
        assert computed_properties["mass_ss_da"] == pytest.approx(bio_mass_ss - 79.98, rel=1e-3)

    def test_short_dna_gc_content(self):
        computed_properties = compute_dna_properties(SHORT_DNA)
        bio_gc = gc_fraction(SHORT_DNA) * 100
        assert computed_properties["gc_prop"] == pytest.approx(bio_gc, rel=1e-4)

    def test_short_dna_melting_temp(self):
        computed_properties = compute_dna_properties(SHORT_DNA)
        gc = gc_fraction(SHORT_DNA) * 100
        bio_tm = 81.5 + 16.6 * math.log10(0.05) + 0.41 * gc - 675 / len(SHORT_DNA)
        assert computed_properties["tm"] == pytest.approx(bio_tm, rel=1e-4)

    def test_long_dna_mass(self):
        computed_properties = compute_dna_properties(LONG_DNA)
        bio_mass_ds = molecular_weight(
            LONG_DNA, "DNA", circular=False, double_stranded=True
        )
        bio_mass_ss = molecular_weight(
            LONG_DNA, "DNA", circular=False, double_stranded=False
        )
        assert computed_properties["mass_ds_da"] == pytest.approx(bio_mass_ds - 159.96, rel=1e-3)
        assert computed_properties["mass_ss_da"] == pytest.approx(bio_mass_ss - 79.98, rel=1e-3)

    def test_long_dna_gc_content(self):
        computed_properties = compute_dna_properties(LONG_DNA)
        bio_gc = gc_fraction(LONG_DNA) * 100
        assert computed_properties["gc_prop"] == pytest.approx(bio_gc, rel=1e-4)

    def test_long_dna_melting_temp(self):
        computed_properties = compute_dna_properties(LONG_DNA)
        gc = gc_fraction(LONG_DNA) * 100
        bio_tm = 81.5 + 16.6 * math.log10(0.05) + 0.41 * gc - 675 / len(LONG_DNA)
        assert computed_properties["tm"] == pytest.approx(bio_tm, rel=1e-4)


class TestRNAProperties:

    def test_short_rna_mass(self):
        computed_properties = compute_rna_properties(SHORT_RNA)
        bio_mass = molecular_weight(
            SHORT_RNA, "RNA", circular=False, double_stranded=False
        )
        assert computed_properties["mass_da"] == pytest.approx(bio_mass - 79.98, rel=1e-3)

    def test_long_rna_mass(self):
        computed_properties = compute_rna_properties(LONG_RNA)
        bio_mass = molecular_weight(
            LONG_RNA, "RNA", circular=False, double_stranded=False
        )
        assert computed_properties["mass_da"] == pytest.approx(bio_mass - 79.98, rel=1e-3)


class TestProteinProperties:

    def test_short_protein_mass(self):
        computed_properties = compute_protein_properties(SHORT_PROT)
        bio_mass = molecular_weight(SHORT_PROT, "protein", circular=False) / 1000.0
        assert computed_properties["mass_kda"] == pytest.approx(bio_mass, rel=1e-3)

    def test_short_protein_pi(self):
        computed_properties = compute_protein_properties(SHORT_PROT)
        bio_pi = ProteinAnalysis(SHORT_PROT).isoelectric_point()
        assert computed_properties["pi"] == pytest.approx(bio_pi, abs=0.9)

    def test_short_protein_extinction_coefficient(self):
        computed_properties = compute_protein_properties(SHORT_PROT)
        bio_ext_reduced, bio_ext_oxidized = ProteinAnalysis(SHORT_PROT).molar_extinction_coefficient()

        assert computed_properties["ext_coeff_reduced"] == pytest.approx(bio_ext_reduced, rel=1e-3)
        assert computed_properties["ext_coeff_oxidized"] == pytest.approx(bio_ext_oxidized, rel=1e-3)

    def test_long_protein_mass(self):
        computed_properties = compute_protein_properties(LONG_PROT)
        bio_mass = molecular_weight(LONG_PROT, "protein", circular=False) / 1000.0
        assert computed_properties["mass_kda"] == pytest.approx(bio_mass, rel=1e-3)

    def test_long_protein_pi(self):
        computed_properties = compute_protein_properties(LONG_PROT)
        bio_pi = ProteinAnalysis(LONG_PROT).isoelectric_point()

        assert computed_properties["pi"] == pytest.approx(bio_pi, abs=0.4)

    def test_long_protein_extinction_coefficient(self):
        computed_properties = compute_protein_properties(LONG_PROT)
        bio_ext_reduced, bio_ext_oxidized = ProteinAnalysis(LONG_PROT).molar_extinction_coefficient()

        assert computed_properties["ext_coeff_reduced"] == pytest.approx(bio_ext_reduced, rel=1e-3)
        assert computed_properties["ext_coeff_oxidized"] == pytest.approx(bio_ext_oxidized, rel=1e-3)
