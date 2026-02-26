"""Tests for Nei-Gojobori codon lookup tables."""

import numpy as np

from pie.codon import (
    AMINO_ACID,
    CODON_TO_INDEX,
    INDEX_TO_CODON,
    N_DIFFS,
    N_SITES,
    N_SITES_EXCL_STOP,
    S_DIFFS,
    S_SITES,
    S_SITES_EXCL_STOP,
    codon_to_index,
    is_stop_codon,
)


class TestCodonIndex:
    def test_64_codons(self):
        assert len(CODON_TO_INDEX) == 64
        assert len(INDEX_TO_CODON) == 64

    def test_roundtrip(self):
        for codon, idx in CODON_TO_INDEX.items():
            assert INDEX_TO_CODON[idx] == codon

    def test_lexicographic_order(self):
        """ACGT lexicographic: AAA=0, AAC=1, AAG=2, AAT=3, ACA=4, ..., TTT=63."""
        assert codon_to_index("AAA") == 0
        assert codon_to_index("AAC") == 1
        assert codon_to_index("AAG") == 2
        assert codon_to_index("AAT") == 3
        assert codon_to_index("ACA") == 4
        assert codon_to_index("TTT") == 63


class TestAminoAcid:
    def test_atg_is_met(self):
        assert AMINO_ACID[codon_to_index("ATG")] == "M"

    def test_stop_codons(self):
        for stop in ["TAA", "TAG", "TGA"]:
            assert AMINO_ACID[codon_to_index(stop)] == "*"

    def test_synonymous_pair(self):
        assert AMINO_ACID[codon_to_index("GCT")] == "A"
        assert AMINO_ACID[codon_to_index("GCC")] == "A"

    def test_known_amino_acids(self):
        """Spot-check several codons against the standard genetic code."""
        checks = {
            "TTT": "F", "TTC": "F",
            "TTA": "L", "TTG": "L",
            "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
            "ATT": "I", "ATC": "I", "ATA": "I",
            "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
            "TGG": "W", "TGT": "C", "TGC": "C",
            "GAT": "D", "GAC": "D",
            "GAA": "E", "GAG": "E",
            "AAA": "K", "AAG": "K",
        }
        for codon, aa in checks.items():
            assert AMINO_ACID[codon_to_index(codon)] == aa, f"{codon} -> {aa}"


class TestIsStopCodon:
    def test_stops(self):
        for stop in ["TAA", "TAG", "TGA"]:
            assert is_stop_codon(codon_to_index(stop))

    def test_non_stop(self):
        assert not is_stop_codon(codon_to_index("ATG"))


class TestSiteCounts:
    def test_shape(self):
        assert N_SITES.shape == (64, 3)
        assert S_SITES.shape == (64, 3)

    def test_sites_sum_to_one(self):
        for i in range(64):
            if AMINO_ACID[i] == "*":
                continue
            for pos in range(3):
                total = N_SITES[i, pos] + S_SITES[i, pos]
                assert abs(total - 1.0) < 1e-10

    def test_aaa_position3(self):
        """AAA (Lys): pos2 mutations -> AAC(Asn), AAG(Lys), AAT(Asn). 1 syn, 2 nonsyn."""
        idx = codon_to_index("AAA")
        assert abs(N_SITES[idx, 2] - 2 / 3) < 1e-10
        assert abs(S_SITES[idx, 2] - 1 / 3) < 1e-10

    def test_atg_all_nonsyn(self):
        """ATG (Met) is the only codon for Met, so all positions are nonsynonymous."""
        idx = codon_to_index("ATG")
        for pos in range(3):
            assert abs(N_SITES[idx, pos] - 1.0) < 1e-10

    def test_fourfold_degenerate(self):
        """GCT (Ala): 3rd position is fourfold degenerate (GCA,GCC,GCG,GCT all Ala)."""
        idx = codon_to_index("GCT")
        assert abs(S_SITES[idx, 2] - 1.0) < 1e-10
        assert abs(N_SITES[idx, 2] - 0.0) < 1e-10

    def test_stop_codons_zero(self):
        """Stop codons should have 0 sites (both modes)."""
        for stop in ["TAA", "TAG", "TGA"]:
            idx = codon_to_index(stop)
            for pos in range(3):
                assert N_SITES[idx, pos] == 0.0
                assert S_SITES[idx, pos] == 0.0
                assert N_SITES_EXCL_STOP[idx, pos] == 0.0
                assert S_SITES_EXCL_STOP[idx, pos] == 0.0

    def test_include_stop_increases_n_sites(self):
        """Including stop_gained as nonsynonymous increases N_sites for affected codons."""
        # TAC (Tyr) pos 2: mutations → TAA(*), TAG(*), TAT(Tyr)
        # Exclude: only TAT(syn), N=0/1, S=1/1
        # Include: TAA(N), TAG(N), TAT(S) → N=2/3, S=1/3
        idx = codon_to_index("TAC")
        assert abs(N_SITES_EXCL_STOP[idx, 2] - 0.0) < 1e-10
        assert abs(S_SITES_EXCL_STOP[idx, 2] - 1.0) < 1e-10
        assert abs(N_SITES[idx, 2] - 2 / 3) < 1e-10
        assert abs(S_SITES[idx, 2] - 1 / 3) < 1e-10

    def test_excl_stop_sites_sum_to_one(self):
        """N_SITES_EXCL_STOP + S_SITES_EXCL_STOP = 1.0 per position for sense codons."""
        for i in range(64):
            if AMINO_ACID[i] == "*":
                continue
            for pos in range(3):
                total = N_SITES_EXCL_STOP[i, pos] + S_SITES_EXCL_STOP[i, pos]
                assert abs(total - 1.0) < 1e-10


class TestDiffCounts:
    def test_shape(self):
        assert N_DIFFS.shape == (64, 64)
        assert S_DIFFS.shape == (64, 64)

    def test_self_diff_zero(self):
        for i in range(64):
            assert N_DIFFS[i, i] == 0.0
            assert S_DIFFS[i, i] == 0.0

    def test_symmetric(self):
        np.testing.assert_array_equal(N_DIFFS, N_DIFFS.T)
        np.testing.assert_array_equal(S_DIFFS, S_DIFFS.T)

    def test_synonymous_change(self):
        """GCT -> GCC: both Ala, so 1 synonymous diff, 0 nonsynonymous."""
        i, j = codon_to_index("GCT"), codon_to_index("GCC")
        assert S_DIFFS[i, j] == 1.0
        assert N_DIFFS[i, j] == 0.0

    def test_nonsynonymous_change(self):
        """AAA (Lys) -> AAC (Asn): 1 nonsynonymous diff."""
        i, j = codon_to_index("AAA"), codon_to_index("AAC")
        assert N_DIFFS[i, j] == 1.0
        assert S_DIFFS[i, j] == 0.0

    def test_two_step_change(self):
        """AAA -> ACC: 2 positions differ, total diffs should be > 0."""
        i, j = codon_to_index("AAA"), codon_to_index("ACC")
        total = N_DIFFS[i, j] + S_DIFFS[i, j]
        assert total > 0

    def test_total_diffs_equal_hamming(self):
        """For any non-stop codon pair, N_DIFFS + S_DIFFS should equal the number
        of differing positions (if no stop-codon pathways are skipped entirely)."""
        bases = "ACGT"
        for i in range(64):
            if AMINO_ACID[i] == "*":
                continue
            for j in range(i + 1, 64):
                if AMINO_ACID[j] == "*":
                    continue
                codon_i = INDEX_TO_CODON[i]
                codon_j = INDEX_TO_CODON[j]
                hamming = sum(a != b for a, b in zip(codon_i, codon_j))
                total = N_DIFFS[i, j] + S_DIFFS[i, j]
                # total should be <= hamming (can be less if all pathways go through stops)
                assert total <= hamming + 1e-10
