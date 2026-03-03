"""Tests for the core piN/piS diversity engine."""

import numpy as np
from pie.codon import codon_to_index, N_SITES, S_SITES, CODON_TO_INDEX
from pie.diversity import (
    _monomorphic_codon_index,
    build_allele_freq_array,
    compute_codon_diversity,
    compute_gene_diversity,
    CodonResult,
    GeneResult,
)
from pie.vcf import FetchResult, FilterStats, Variant


class TestBuildAlleleFreqArray:
    """Tests for building the (n_codons, 3, 4) allele frequency array."""

    def test_monomorphic(self):
        """All-reference codons: each position has freq 1.0 for ref allele."""
        codons = ["ATG", "AAA", "GCT"]
        positions = [("chr1", 0, 1, 2), ("chr1", 3, 4, 5), ("chr1", 6, 7, 8)]
        freqs = build_allele_freq_array(codons, positions, [])
        assert freqs.shape == (3, 3, 4)
        # ATG: A=0 at pos0, T=3 at pos1, G=2 at pos2
        assert freqs[0, 0, 0] == 1.0  # A at pos 0
        assert freqs[0, 1, 3] == 1.0  # T at pos 1
        assert freqs[0, 2, 2] == 1.0  # G at pos 2
        # AAA: A=0 at all positions
        assert freqs[1, 0, 0] == 1.0
        assert freqs[1, 1, 0] == 1.0
        assert freqs[1, 2, 0] == 1.0
        # GCT: G=2, C=1, T=3
        assert freqs[2, 0, 2] == 1.0
        assert freqs[2, 1, 1] == 1.0
        assert freqs[2, 2, 3] == 1.0

    def test_one_variant(self):
        """Single variant reduces ref freq, adds alt freq."""
        codons = ["AAA"]
        positions = [("chr1", 0, 1, 2)]
        variants = [Variant(pos=2, ref="A", alt="G", freq=0.1, depth=100)]
        freqs = build_allele_freq_array(codons, positions, variants)
        assert abs(freqs[0, 2, 0] - 0.9) < 1e-10  # A ref reduced
        assert abs(freqs[0, 2, 2] - 0.1) < 1e-10  # G alt added

    def test_two_variants_same_position(self):
        """Two biallelic records at same position -> 3 alleles."""
        codons = ["AAA"]
        positions = [("chr1", 0, 1, 2)]
        variants = [
            Variant(pos=2, ref="A", alt="G", freq=0.1, depth=100),
            Variant(pos=2, ref="A", alt="C", freq=0.05, depth=100),
        ]
        freqs = build_allele_freq_array(codons, positions, variants)
        assert abs(freqs[0, 2, 0] - 0.85) < 1e-10  # A (1.0 - 0.1 - 0.05)
        assert abs(freqs[0, 2, 2] - 0.1) < 1e-10   # G
        assert abs(freqs[0, 2, 1] - 0.05) < 1e-10   # C

    def test_ambiguous_base_n_defensive_guard(self):
        """Non-ACGT base (N) in codon: break leaves remaining positions zeroed (Issue #2).

        Upstream (compute_gene_diversity) filters codons with N before calling
        this function. The defensive guard here breaks on N, leaving the N
        position and any subsequent positions as zeros.
        """
        codons = ["ATN", "AAA"]
        positions = [("chr1", 0, 1, 2), ("chr1", 3, 4, 5)]
        freqs = build_allele_freq_array(codons, positions, [])
        # First codon: A and T are set, N triggers break -> pos 2 all zeros
        assert freqs[0, 0, 0] == 1.0  # A set before break
        assert freqs[0, 1, 3] == 1.0  # T set before break
        assert np.all(freqs[0, 2] == 0.0)  # N position zeroed
        # Second codon is clean
        assert freqs[1, 0, 0] == 1.0  # A
        assert freqs[1, 1, 0] == 1.0  # A
        assert freqs[1, 2, 0] == 1.0  # A

    def test_minus_strand_complement(self):
        """On - strand, variant REF/ALT must be complemented."""
        # CDS-sense codon is "ATG", genomic positions are reversed
        # At genomic pos 10, CDS-sense base is 'G' (3rd codon pos)
        # On - strand, genomic base = complement('G') = 'C'
        # VCF variant: REF=C, ALT=A at pos 10
        # In CDS-sense: complement(C)=G (ref), complement(A)=T (alt)
        codons = ["ATG"]
        positions = [("chr1", 12, 11, 10)]  # reversed for - strand
        variants = [Variant(pos=10, ref="C", alt="A", freq=0.3, depth=100)]
        freqs = build_allele_freq_array(codons, positions, variants, strand="-")
        # CDS-sense: ref=G (idx=2), alt=T (idx=3) at codon position 2
        assert abs(freqs[0, 2, 2] - 0.7) < 1e-10   # G (ref) reduced
        assert abs(freqs[0, 2, 3] - 0.3) < 1e-10   # T (alt) added


class TestComputeCodonDiversity:
    """Tests for single-codon diversity computation."""

    def test_monomorphic(self):
        """Monomorphic codon: no diffs, only site counts."""
        freq = np.zeros((3, 4))
        freq[0, 0] = 1.0  # A
        freq[1, 0] = 1.0  # A
        freq[2, 0] = 1.0  # A -> AAA
        result = compute_codon_diversity(freq)
        assert result["N_diffs"] == 0.0
        assert result["S_diffs"] == 0.0
        idx = codon_to_index("AAA")
        assert abs(result["N_sites"] - float(N_SITES[idx].sum())) < 1e-10
        assert abs(result["S_sites"] - float(S_SITES[idx].sum())) < 1e-10

    def test_synonymous_variant(self):
        """GCT with T->A at pos3 (GCT->GCA, both Ala): S_diffs > 0, N_diffs = 0."""
        freq = np.zeros((3, 4))
        freq[0, 2] = 1.0  # G
        freq[1, 1] = 1.0  # C
        freq[2, 3] = 0.8  # T at 0.8
        freq[2, 0] = 0.2  # A at 0.2
        result = compute_codon_diversity(freq)
        assert result["S_diffs"] > 0
        assert result["N_diffs"] == 0.0
        # Expected: S_diffs = 2 * 0.8 * 0.2 * 1.0 = 0.32
        # (GCT->GCA is 1 synonymous diff)
        assert abs(result["S_diffs"] - 0.32) < 1e-10

    def test_nonsynonymous_variant(self):
        """AAA with A->C at pos3 (AAA->AAC, Lys->Asn): N_diffs > 0, S_diffs = 0."""
        freq = np.zeros((3, 4))
        freq[0, 0] = 1.0  # A
        freq[1, 0] = 1.0  # A
        freq[2, 0] = 0.9  # A
        freq[2, 1] = 0.1  # C
        result = compute_codon_diversity(freq)
        assert result["N_diffs"] > 0
        assert result["S_diffs"] == 0.0
        assert abs(result["N_diffs"] - 0.18) < 1e-10  # 2 * 0.9 * 0.1 * 1.0

    def test_stop_gained_as_nonsynonymous(self):
        """TGG->TGA (stop_gained): counts as nonsynonymous by default."""
        freq = np.zeros((3, 4))
        freq[0, 3] = 1.0  # T
        freq[1, 2] = 1.0  # G
        freq[2, 2] = 0.95  # G (TGG=Trp)
        freq[2, 0] = 0.05  # A (TGA=Stop)
        result = compute_codon_diversity(freq)
        # TGG->TGA is nonsynonymous (Trp->Stop)
        # N_diffs = 2 * 0.95 * 0.05 * 1.0 = 0.095
        assert abs(result["N_diffs"] - 0.095) < 1e-10
        assert result["S_diffs"] == 0.0

    def test_stop_codon_excluded_legacy(self):
        """TGG->TGA with exclude_stops=True: stop removed, old behavior."""
        freq = np.zeros((3, 4))
        freq[0, 3] = 1.0  # T
        freq[1, 2] = 1.0  # G
        freq[2, 2] = 0.95  # G (TGG=Trp)
        freq[2, 0] = 0.05  # A (TGA=Stop)
        result = compute_codon_diversity(freq, exclude_stops=True)
        # After stop removal and renormalization, only TGG remains
        assert result["N_diffs"] == 0.0
        assert result["S_diffs"] == 0.0

    def test_weighted_site_counts(self):
        """Site counts should be weighted average of constituent codons."""
        # GCT(Ala) freq=0.8, GCC(Ala) freq=0.2 -> weighted N_sites and S_sites
        freq = np.zeros((3, 4))
        freq[0, 2] = 1.0  # G
        freq[1, 1] = 1.0  # C
        freq[2, 3] = 0.8  # T (GCT)
        freq[2, 1] = 0.2  # C (GCC)
        result = compute_codon_diversity(freq)
        idx_gct = codon_to_index("GCT")
        idx_gcc = codon_to_index("GCC")
        expected_n = 0.8 * float(N_SITES[idx_gct].sum()) + 0.2 * float(N_SITES[idx_gcc].sum())
        expected_s = 0.8 * float(S_SITES[idx_gct].sum()) + 0.2 * float(S_SITES[idx_gcc].sum())
        assert abs(result["N_sites"] - expected_n) < 1e-10
        assert abs(result["S_sites"] - expected_s) < 1e-10


class TestComputeGeneDiversity:
    """Integration tests with the full test dataset."""

    def test_gene1(self, ref_fasta, gff3_file, vcf_file):
        """Gene1: known piN and piS from hand calculations."""
        from pie.reference import ReferenceGenome
        from pie.annotation import parse_annotations
        from pie.vcf import VariantReader

        genes = parse_annotations(gff3_file)
        gene1 = [g for g in genes if "gene1" in g.gene_id.lower()][0]

        with ReferenceGenome(ref_fasta) as ref, \
             VariantReader(vcf_file, min_freq=0.0, min_depth=0, min_qual=0) as vcf:
            result = compute_gene_diversity(gene1, ref, vcf)

        assert isinstance(result, GeneResult)
        assert result.n_codons == 29  # 30 total - 1 stop
        assert result.n_poly_codons == 2
        assert result.n_variants == 2

        # Hand-calculated expected values:
        # N_sites = ATG(3.0) + GAT(8/3) + 27*GCT(2.0) = 59.6667
        assert abs(result.N_sites - 59.6667) < 0.001
        # S_sites = ATG(0.0) + GAT(1/3) + 27*GCT(1.0) = 27.3333
        assert abs(result.S_sites - 27.3333) < 0.001
        # N_diffs = 2 * 0.70 * 0.30 = 0.42 (variant at pos 7, G->A, GAT->AAT)
        assert abs(result.N_diffs - 0.42) < 1e-10
        # S_diffs = 2 * 0.80 * 0.20 = 0.32 (variant at pos 6, T->C, GCT->GCC)
        assert abs(result.S_diffs - 0.32) < 1e-10
        # piN = 0.42 / 59.6667 ~ 0.007039
        assert abs(result.piN - 0.007039) < 0.0001
        # piS = 0.32 / 27.3333 ~ 0.011707
        assert abs(result.piS - 0.011707) < 0.0001

    def test_gene2(self, ref_fasta, gff3_file, vcf_file):
        """Gene2: multi-exon, 1 nonsynonymous variant."""
        from pie.reference import ReferenceGenome
        from pie.annotation import parse_annotations
        from pie.vcf import VariantReader

        genes = parse_annotations(gff3_file)
        gene2 = [g for g in genes if "gene2" in g.gene_id.lower()][0]

        with ReferenceGenome(ref_fasta) as ref, \
             VariantReader(vcf_file, min_freq=0.0, min_depth=0, min_qual=0) as vcf:
            result = compute_gene_diversity(gene2, ref, vcf)

        assert result.n_codons == 32  # 33 - 1 stop
        assert result.n_poly_codons == 1
        # N_diffs = 2 * 0.60 * 0.40 = 0.48
        assert abs(result.N_diffs - 0.48) < 1e-10
        assert result.S_diffs == 0.0
        assert result.piS == 0.0
        assert result.piN_piS is None  # piS = 0

    def test_gene3_with_qual_filter(self, ref_fasta, gff3_file, vcf_file):
        """Gene3: minus strand, variant filtered by QUAL=15 < 20."""
        from pie.reference import ReferenceGenome
        from pie.annotation import parse_annotations
        from pie.vcf import VariantReader

        genes = parse_annotations(gff3_file)
        gene3 = [g for g in genes if "gene3" in g.gene_id.lower()][0]

        with ReferenceGenome(ref_fasta) as ref, \
             VariantReader(vcf_file, min_freq=0.0, min_depth=0, min_qual=20) as vcf:
            result = compute_gene_diversity(gene3, ref, vcf)

        assert result.n_variants == 0
        assert result.N_diffs == 0.0
        assert result.S_diffs == 0.0
        assert result.piN == 0.0
        assert result.piS == 0.0

    def test_gene3_without_qual_filter(self, ref_fasta, gff3_file, vcf_file):
        """Gene3: minus strand, variant passes when QUAL filter is disabled."""
        from pie.reference import ReferenceGenome
        from pie.annotation import parse_annotations
        from pie.vcf import VariantReader

        genes = parse_annotations(gff3_file)
        gene3 = [g for g in genes if "gene3" in g.gene_id.lower()][0]

        with ReferenceGenome(ref_fasta) as ref, \
             VariantReader(vcf_file, min_freq=0.0, min_depth=0, min_qual=0) as vcf:
            result = compute_gene_diversity(gene3, ref, vcf)

        assert result.n_variants == 1
        assert result.n_poly_codons == 1
        # Synonym variant: S_diffs = 2 * 0.5 * 0.5 = 0.50
        assert abs(result.S_diffs - 0.50) < 1e-10
        assert result.N_diffs == 0.0

    def test_gene_result_properties(self, ref_fasta, gff3_file, vcf_file):
        """Verify GeneResult properties piN, piS, piN_piS."""
        from pie.reference import ReferenceGenome
        from pie.annotation import parse_annotations
        from pie.vcf import VariantReader

        genes = parse_annotations(gff3_file)
        gene1 = [g for g in genes if "gene1" in g.gene_id.lower()][0]

        with ReferenceGenome(ref_fasta) as ref, \
             VariantReader(vcf_file, min_freq=0.0, min_depth=0, min_qual=0) as vcf:
            result = compute_gene_diversity(gene1, ref, vcf)

        # piN_piS should be piN / piS
        assert result.piN_piS is not None
        expected_ratio = result.piN / result.piS
        assert abs(result.piN_piS - expected_ratio) < 1e-10

    def test_codon_results_populated(self, ref_fasta, gff3_file, vcf_file):
        """Verify codon-level results are stored correctly."""
        from pie.reference import ReferenceGenome
        from pie.annotation import parse_annotations
        from pie.vcf import VariantReader

        genes = parse_annotations(gff3_file)
        gene1 = [g for g in genes if "gene1" in g.gene_id.lower()][0]

        with ReferenceGenome(ref_fasta) as ref, \
             VariantReader(vcf_file, min_freq=0.0, min_depth=0, min_qual=0) as vcf:
            result = compute_gene_diversity(gene1, ref, vcf)

        assert len(result.codon_results) == 29  # excluding stop
        # All codon results should have valid site counts
        for cr in result.codon_results:
            assert cr.N_sites >= 0
            assert cr.S_sites >= 0
            assert isinstance(cr, CodonResult)


class TestMonomorphicCodonIndex:
    """Tests for _monomorphic_codon_index helper."""

    def test_reference_codon(self):
        """Reference-only codon returns the correct index."""
        freq = np.zeros((3, 4))
        freq[0, 0] = 1.0  # A
        freq[1, 0] = 1.0  # A
        freq[2, 0] = 1.0  # A -> AAA
        assert _monomorphic_codon_index(freq) == CODON_TO_INDEX["AAA"]

    def test_fixed_alt_codon(self):
        """Fixed alternate allele returns the derived codon index."""
        freq = np.zeros((3, 4))
        freq[0, 3] = 1.0  # T (was A in ref)
        freq[1, 0] = 1.0  # A
        freq[2, 0] = 1.0  # A -> TAA
        assert _monomorphic_codon_index(freq) == CODON_TO_INDEX["TAA"]

    def test_all_zeros_returns_none(self):
        """All-zero position returns None."""
        freq = np.zeros((3, 4))
        freq[0, 0] = 1.0
        freq[1, 0] = 1.0
        # pos 2 all zeros
        assert _monomorphic_codon_index(freq) is None

    def test_polymorphic_returns_none(self):
        """Polymorphic position (2 alleles) returns None."""
        freq = np.zeros((3, 4))
        freq[0, 0] = 0.8
        freq[0, 3] = 0.2  # two alleles at pos 0
        freq[1, 0] = 1.0
        freq[2, 0] = 1.0
        assert _monomorphic_codon_index(freq) is None


class TestFixedAltMonomorphicSiteCounts:
    """Regression tests for issue #4: fixed non-reference monomorphic sites
    must use the derived codon for N_sites/S_sites, not the reference codon."""

    def test_fixed_alt_uses_derived_site_counts(self):
        """A codon fixed for ALT should get site counts from the derived codon,
        not the reference codon."""
        # Reference: AAA (Lys, N_sites=2.6667)
        # Fixed derived: ACA (Thr, N_sites=2.0000) — pos 1 A->C
        codons = ["AAA"]  # reference codon from FASTA
        positions = [("chr1", 0, 1, 2)]
        variants = [Variant(pos=1, ref="A", alt="C", freq=1.0, depth=100)]
        freq_array = build_allele_freq_array(codons, positions, variants)

        # Verify freq_array encodes ACA (C at pos 1)
        assert freq_array[0, 1, 1] == 1.0  # C
        assert freq_array[0, 1, 0] == 0.0  # A gone

        # Monomorphic: poly_mask should be False
        poly_mask = (freq_array > 0).sum(axis=2).max(axis=1) > 1
        assert not poly_mask[0]

        # The actual codon index should be ACA, not AAA
        actual_idx = _monomorphic_codon_index(freq_array[0])
        assert actual_idx == CODON_TO_INDEX["ACA"]
        assert actual_idx != CODON_TO_INDEX["AAA"]

        # Site counts must differ between AAA and ACA
        ref_n = float(N_SITES[CODON_TO_INDEX["AAA"]].sum())
        derived_n = float(N_SITES[CODON_TO_INDEX["ACA"]].sum())
        assert abs(ref_n - derived_n) > 0.1  # delta is 0.6667

    def test_reference_monomorphic_unchanged(self):
        """Reference-only codons still use reference site counts (no regression)."""
        codons = ["GCT"]
        positions = [("chr1", 0, 1, 2)]
        freq_array = build_allele_freq_array(codons, positions, [])

        actual_idx = _monomorphic_codon_index(freq_array[0])
        assert actual_idx == CODON_TO_INDEX["GCT"]

    def test_gene_with_fixed_alt_variant(self):
        """Integration: compute_gene_diversity uses derived codon site counts
        for a gene with a fixed non-reference allele."""
        from unittest.mock import MagicMock

        from pie.annotation import GeneModel

        # Reference: AAA GCT TAA (3 codons, last is stop)
        # Variant: pos 1 A->C at freq=1.0 → AAA becomes ACA (Thr)
        # AAA N_sites=2.6667 vs ACA N_sites=2.0000 (delta=0.6667)
        gene = GeneModel(
            gene_id="test_fixed", transcript_id="tx1",
            chrom="chr1", start=0, end=9, strand="+",
            cds_exons=[("chr1", 0, 9)],
        )

        ref = MagicMock()
        ref.extract_codons.return_value = ["AAA", "GCT", "TAA"]
        ref.codon_genomic_positions.return_value = [
            ("chr1", 0, 1, 2), ("chr1", 3, 4, 5), ("chr1", 6, 7, 8),
        ]

        # pos 1: A->C at freq=1.0 (AAA -> ACA, fixed derived)
        vcf = MagicMock()
        vcf.fetch.return_value = FetchResult(
            [Variant(pos=1, ref="A", alt="C", freq=1.0, depth=100)],
            FilterStats(),
        )

        result = compute_gene_diversity(gene, ref, vcf)

        assert result.n_codons == 2  # AAA(->ACA) + GCT, TAA skipped
        assert result.n_poly_codons == 0
        assert result.N_diffs == 0.0
        assert result.S_diffs == 0.0

        # First codon should have ACA site counts, not AAA
        cr0 = result.codon_results[0]
        expected_n = float(N_SITES[CODON_TO_INDEX["ACA"]].sum())
        expected_s = float(S_SITES[CODON_TO_INDEX["ACA"]].sum())
        assert abs(cr0.N_sites - expected_n) < 1e-10
        assert abs(cr0.S_sites - expected_s) < 1e-10

        # Confirm it's NOT the reference AAA site counts
        ref_n = float(N_SITES[CODON_TO_INDEX["AAA"]].sum())
        assert abs(cr0.N_sites - ref_n) > 0.1  # delta is 0.6667
