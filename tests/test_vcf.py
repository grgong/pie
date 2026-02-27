import os
import shutil
from pie.vcf import ensure_indexed, VariantReader, IndividualVariantReader


class TestEnsureIndexed:
    def test_already_indexed(self, vcf_file):
        result = ensure_indexed(vcf_file)
        assert result.endswith(".vcf.gz")

    def test_plain_vcf_gets_indexed(self, plain_vcf_file, tmp_path):
        tmp_vcf = str(tmp_path / "test.vcf")
        shutil.copy(plain_vcf_file, tmp_vcf)
        result = ensure_indexed(tmp_vcf)
        assert result.endswith(".vcf.gz")
        assert os.path.exists(result + ".tbi")


class TestVariantReader:
    def test_context_manager(self, vcf_file):
        with VariantReader(vcf_file, min_freq=0.0, min_depth=0, min_qual=0) as reader:
            assert reader is not None

    def test_fetch_gene1_variants(self, vcf_file):
        """Gene1 region has 2 variants: pos 5 (T->C) and pos 6 (G->A) (0-based)."""
        with VariantReader(vcf_file, min_freq=0.0, min_depth=0, min_qual=0) as reader:
            variants = reader.fetch("chr1", 0, 90)
            assert len(variants) == 2
            # Check first variant
            v1 = variants[0]
            assert v1.pos == 5  # VCF pos 6 (1-based) -> 5 (0-based)
            assert v1.ref == "T"
            assert v1.alt == "C"
            assert abs(v1.freq - 0.20) < 1e-6
            assert v1.depth == 100

    def test_fetch_gene2_variant(self, vcf_file):
        with VariantReader(vcf_file, min_freq=0.0, min_depth=0, min_qual=0) as reader:
            variants = reader.fetch("chr1", 100, 220)
            assert len(variants) == 1
            assert variants[0].pos == 194  # VCF 195 -> 0-based 194
            assert abs(variants[0].freq - 0.40) < 1e-6

    def test_min_freq_filter(self, vcf_file):
        with VariantReader(vcf_file, min_freq=0.99, min_depth=0, min_qual=0) as reader:
            variants = reader.fetch("chr1", 0, 350)
            assert len(variants) == 0

    def test_min_depth_filter(self, vcf_file):
        with VariantReader(vcf_file, min_freq=0.0, min_depth=999999, min_qual=0) as reader:
            variants = reader.fetch("chr1", 0, 350)
            assert len(variants) == 0

    def test_qual_filter(self, vcf_file):
        """Variant 4 at pos 297 has QUAL=15, should be filtered at min_qual=20."""
        with VariantReader(vcf_file, min_freq=0.0, min_depth=0, min_qual=20) as reader:
            variants = reader.fetch("chr1", 230, 311)
            assert len(variants) == 0  # variant 4 filtered

    def test_qual_filter_disabled(self, vcf_file):
        """With min_qual=0, variant 4 should pass."""
        with VariantReader(vcf_file, min_freq=0.0, min_depth=0, min_qual=0) as reader:
            variants = reader.fetch("chr1", 230, 311)
            assert len(variants) == 1
            assert variants[0].pos == 296


class TestMultiallelicFiltering:
    def test_default_skips_multiallelic(self, multiallelic_vcf_file):
        """Default behavior: positions with >1 ALT allele are discarded."""
        with VariantReader(multiallelic_vcf_file, min_freq=0.0, min_depth=0,
                           min_qual=0) as reader:
            variants = reader.fetch("chr1", 0, 350)
            # pos 7 has 2 ALT alleles (decomposed) -> skipped
            # Only pos 6 and pos 195 remain
            assert len(variants) == 2
            positions = [v.pos for v in variants]
            assert 5 in positions   # pos 6 (1-based) -> 5 (0-based)
            assert 194 in positions  # pos 195 -> 194
            assert 6 not in positions  # pos 7 -> 6, should be skipped

    def test_keep_multiallelic_preserves_all(self, multiallelic_vcf_file):
        """With keep_multiallelic=True, multiallelic sites are merged and kept."""
        with VariantReader(multiallelic_vcf_file, min_freq=0.0, min_depth=0,
                           min_qual=0, keep_multiallelic=True) as reader:
            variants = reader.fetch("chr1", 0, 350)
            # pos 7 has 2 ALT alleles -> merged and kept
            assert len(variants) == 4
            positions = [v.pos for v in variants]
            assert positions.count(6) == 2  # two ALT alleles at pos 7 (0-based 6)

    def test_multiallelic_still_skipped_after_freq_filter(self, multiallelic_vcf_file):
        """Multiallelic site must be skipped even when one ALT is filtered by min_freq.

        pos 7 decomposed records: AD=50,30 (freq=0.375) and AD=50,20 (freq=0.286).
        With min_freq=0.30, ALT=C (0.286) is filtered out but the site is still
        multiallelic and must be skipped.
        """
        with VariantReader(multiallelic_vcf_file, min_freq=0.30, min_depth=0,
                           min_qual=0) as reader:
            variants = reader.fetch("chr1", 0, 350)
            positions = [v.pos for v in variants]
            assert 6 not in positions, "multiallelic pos 7 should be skipped even if one ALT filtered"

    def test_keep_multiallelic_freq_uses_complete_depth(self, multiallelic_vcf_file):
        """Frequencies must be computed from ALL allele depths, not just filtered ones.

        pos 7 decomposed: AD=50,30 (A) and AD=50,20 (C).
        total_depth = 50 + 30 + 20 = 100.
        With keep_multiallelic + min_freq=0.25:
        - ALT=A freq = 30/100 = 0.30 -> passes
        - ALT=C freq = 20/100 = 0.20 -> filtered
        BUG (old): if C is pre-filtered, depth=50+30=80, freq(A)=30/80=0.375 (inflated)
        """
        with VariantReader(multiallelic_vcf_file, min_freq=0.25, min_depth=0,
                           min_qual=0, keep_multiallelic=True) as reader:
            variants = reader.fetch("chr1", 0, 350)
            ma_variants = [v for v in variants if v.pos == 6]
            assert len(ma_variants) == 1  # only ALT=A passes min_freq
            v = ma_variants[0]
            assert v.alt == "A"
            assert v.depth == 100  # ref(50) + A(30) + C(20)
            assert abs(v.freq - 0.30) < 1e-6  # 30/100, not 30/80

    def test_keep_multiallelic_depth_uses_complete_depth(self, multiallelic_vcf_file):
        """min_depth must not pre-filter alleles before multiallelic merge.

        pos 7 decomposed: AD=50,30 (A, per-record depth=80) and
        AD=50,20 (C, per-record depth=70).
        With keep_multiallelic + min_depth=75:
        - Record C has per-record depth 70 < 75
        - But merged total_depth = 50 + 30 + 20 = 100 >= 75
        - ALT=A freq = 30/100 = 0.30
        - ALT=C freq = 20/100 = 0.20
        BUG (old): C pre-filtered by depth, merge sees only A:
        total_depth=50+30=80, freq(A)=30/80=0.375 (inflated)
        """
        with VariantReader(multiallelic_vcf_file, min_freq=0.0, min_depth=75,
                           min_qual=0, keep_multiallelic=True) as reader:
            variants = reader.fetch("chr1", 0, 350)
            ma_variants = [v for v in variants if v.pos == 6]
            assert len(ma_variants) == 2  # both ALTs survive
            for v in ma_variants:
                assert v.depth == 100  # ref(50) + A(30) + C(20)
            va = next(v for v in ma_variants if v.alt == "A")
            vc = next(v for v in ma_variants if v.alt == "C")
            assert abs(va.freq - 0.30) < 1e-6  # 30/100, not 30/80
            assert abs(vc.freq - 0.20) < 1e-6  # 20/100

    def test_non_decomposed_multiallelic_skipped(self, multiallelic_inline_vcf_file):
        """Non-decomposed multiallelic (single line ALT=A,C) is also skipped."""
        with VariantReader(multiallelic_inline_vcf_file, min_freq=0.0, min_depth=0,
                           min_qual=0) as reader:
            variants = reader.fetch("chr1", 0, 350)
            assert len(variants) == 2
            positions = [v.pos for v in variants]
            assert 5 in positions
            assert 194 in positions
            assert 6 not in positions


class TestIndividualVariantReader:
    def test_basic_gt_frequency(self, individual_vcf_file):
        """All 4 samples, no filtering. Verify GT-derived frequencies."""
        with IndividualVariantReader(
            individual_vcf_file,
            samples=["S1", "S2", "S3", "S4"],
            min_freq=0.0, min_qual=0.0, pass_only=False,
            keep_multiallelic=False, min_call_rate=0.0, min_an=0,
        ) as reader:
            variants = reader.fetch("chr1", 0, 350)
            # 4 biallelic sites, all pass with no filters
            assert len(variants) == 4

            v6 = next(v for v in variants if v.pos == 5)
            assert v6.ref == "T" and v6.alt == "C"
            assert abs(v6.freq - 2 / 6) < 1e-6   # AC=2, AN=6
            assert v6.depth == 6                    # AN
            assert abs(v6.call_rate - 0.75) < 1e-6  # 3/4

            v7 = next(v for v in variants if v.pos == 6)
            assert abs(v7.freq - 2 / 8) < 1e-6
            assert v7.depth == 8
            assert abs(v7.call_rate - 1.0) < 1e-6

            v195 = next(v for v in variants if v.pos == 194)
            assert abs(v195.freq - 5 / 8) < 1e-6  # 3 hets(AC=3) + 1 hom_alt(AC=2) = 5 alt alleles
            assert v195.depth == 8

            v297 = next(v for v in variants if v.pos == 296)
            assert abs(v297.freq - 1 / 4) < 1e-6
            assert v297.depth == 4
            assert abs(v297.call_rate - 0.50) < 1e-6

    def test_context_manager(self, individual_vcf_file):
        with IndividualVariantReader(
            individual_vcf_file, samples=["S1", "S2", "S3", "S4"],
            min_freq=0.0, min_qual=0.0, pass_only=False,
            keep_multiallelic=False, min_call_rate=0.0, min_an=0,
        ) as reader:
            assert reader is not None
            assert reader.n_samples == 4

    def test_min_call_rate_filter(self, individual_vcf_file):
        """min_call_rate=0.8 filters pos 6 (0.75) and pos 297 (0.50)."""
        with IndividualVariantReader(
            individual_vcf_file, samples=["S1", "S2", "S3", "S4"],
            min_freq=0.0, min_qual=0.0, pass_only=False,
            keep_multiallelic=False, min_call_rate=0.8, min_an=0,
        ) as reader:
            variants = reader.fetch("chr1", 0, 350)
            positions = [v.pos for v in variants]
            assert 5 not in positions   # pos 6: call_rate=0.75 < 0.8
            assert 6 in positions       # pos 7: call_rate=1.0
            assert 194 in positions     # pos 195: call_rate=1.0
            assert 296 not in positions  # pos 297: call_rate=0.50 < 0.8
            assert len(variants) == 2

    def test_min_an_filter(self, individual_vcf_file):
        """min_an=6 filters pos 297 (AN=4) but keeps pos 6 (AN=6)."""
        with IndividualVariantReader(
            individual_vcf_file, samples=["S1", "S2", "S3", "S4"],
            min_freq=0.0, min_qual=0.0, pass_only=False,
            keep_multiallelic=False, min_call_rate=0.0, min_an=6,
        ) as reader:
            variants = reader.fetch("chr1", 0, 350)
            positions = [v.pos for v in variants]
            assert 5 in positions       # AN=6 >= 6
            assert 6 in positions       # AN=8 >= 6
            assert 194 in positions     # AN=8 >= 6
            assert 296 not in positions  # AN=4 < 6

    def test_min_freq_on_gt_derived_af(self, individual_vcf_file):
        """min_freq=0.30 filters pos 7 (freq=0.25) and pos 297 (freq=0.25)."""
        with IndividualVariantReader(
            individual_vcf_file, samples=["S1", "S2", "S3", "S4"],
            min_freq=0.30, min_qual=0.0, pass_only=False,
            keep_multiallelic=False, min_call_rate=0.0, min_an=0,
        ) as reader:
            variants = reader.fetch("chr1", 0, 350)
            positions = [v.pos for v in variants]
            assert 5 in positions       # freq=2/6=0.333 >= 0.30
            assert 6 not in positions   # freq=2/8=0.25 < 0.30
            assert 194 in positions     # freq=5/8=0.625 >= 0.30
            assert 296 not in positions  # freq=1/4=0.25 < 0.30

    def test_qual_filter(self, individual_vcf_file):
        """min_qual=20 filters pos 297 (QUAL=15)."""
        with IndividualVariantReader(
            individual_vcf_file, samples=["S1", "S2", "S3", "S4"],
            min_freq=0.0, min_qual=20.0, pass_only=False,
            keep_multiallelic=False, min_call_rate=0.0, min_an=0,
        ) as reader:
            variants = reader.fetch("chr1", 0, 350)
            positions = [v.pos for v in variants]
            assert 296 not in positions  # QUAL=15 < 20

    def test_sample_subset(self, individual_vcf_file):
        """Selecting only S1,S2: pos 6 T>C has AC=1/AN=4, call_rate=1.0."""
        with IndividualVariantReader(
            individual_vcf_file, samples=["S1", "S2"],
            min_freq=0.0, min_qual=0.0, pass_only=False,
            keep_multiallelic=False, min_call_rate=0.0, min_an=0,
        ) as reader:
            assert reader.n_samples == 2
            variants = reader.fetch("chr1", 0, 100)
            v6 = next(v for v in variants if v.pos == 5)
            # S1: 0/1 (AC=1), S2: 0/0 (AC=0) -> AC=1, AN=4
            assert abs(v6.freq - 1 / 4) < 1e-6
            assert v6.depth == 4
            assert abs(v6.call_rate - 1.0) < 1e-6  # both called

    def test_all_missing_site_skipped(self, individual_vcf_file):
        """If all selected samples are missing, the site is skipped."""
        with IndividualVariantReader(
            individual_vcf_file, samples=["S4"],
            min_freq=0.0, min_qual=0.0, pass_only=False,
            keep_multiallelic=False, min_call_rate=0.0, min_an=0,
        ) as reader:
            variants = reader.fetch("chr1", 0, 10)
            # pos 6: S4 is ./., so 0 called -> skipped
            positions = [v.pos for v in variants]
            assert 5 not in positions
