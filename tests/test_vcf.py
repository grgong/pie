import os
import shutil

import pytest

from pie.vcf import ensure_indexed, get_vcf_contigs, VariantReader, IndividualVariantReader
from tests.helpers import bgzip_and_index


_GT_FORMAT = '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
_AD_FORMAT = '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depth">'


def _write_vcf(tmp_path, name, samples, rows, fmt=_GT_FORMAT):
    """Minimal single-contig VCF with one FORMAT field. Returns .vcf.gz path."""
    p = tmp_path / name
    p.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=1000>\n"
        + fmt + "\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        + "\t".join(samples) + "\n"
        + "".join(row + "\n" for row in rows)
    )
    return bgzip_and_index(p)


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


class TestGetVcfContigs:
    def test_returns_frozenset(self, vcf_file):
        contigs = get_vcf_contigs(vcf_file)
        assert isinstance(contigs, frozenset)
        assert "chr1" in contigs

    def test_mismatch_vcf_contigs(self, mismatch_vcf_file):
        contigs = get_vcf_contigs(mismatch_vcf_file)
        assert "1" in contigs
        assert "chr1" not in contigs


class TestMissingContigTracking:
    def test_fetch_missing_contig_returns_empty(self, vcf_file):
        """Fetching from absent contig returns [] and tracks it (Issue #7)."""
        with VariantReader(vcf_file, min_freq=0.0, min_depth=0, min_qual=0) as reader:
            variants = reader.fetch("nonexistent", 0, 100).variants
            assert variants == []
            assert "nonexistent" in reader._missing_contigs

    def test_individual_reader_missing_contig(self, individual_vcf_file):
        """IndividualVariantReader also tracks missing contigs."""
        with IndividualVariantReader(
            individual_vcf_file, samples=["S1", "S2", "S3", "S4"],
            min_freq=0.0, min_qual=0.0, pass_only=False,
            keep_multiallelic=False, min_call_rate=0.0, min_an=0,
        ) as reader:
            variants = reader.fetch("chrX", 0, 100).variants
            assert variants == []
            assert "chrX" in reader._missing_contigs


class TestVariantReader:
    def test_context_manager(self, vcf_file):
        with VariantReader(vcf_file, min_freq=0.0, min_depth=0, min_qual=0) as reader:
            assert reader is not None

    def test_fetch_gene1_variants(self, vcf_file):
        """Gene1 region has 2 variants: pos 5 (T->C) and pos 6 (G->A) (0-based)."""
        with VariantReader(vcf_file, min_freq=0.0, min_depth=0, min_qual=0) as reader:
            variants = reader.fetch("chr1", 0, 90).variants
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
            variants = reader.fetch("chr1", 100, 220).variants
            assert len(variants) == 1
            assert variants[0].pos == 194  # VCF 195 -> 0-based 194
            assert abs(variants[0].freq - 0.40) < 1e-6

    def test_min_freq_filter(self, vcf_file):
        with VariantReader(vcf_file, min_freq=0.99, min_depth=0, min_qual=0) as reader:
            variants = reader.fetch("chr1", 0, 350).variants
            assert len(variants) == 0

    def test_min_depth_filter(self, vcf_file):
        with VariantReader(vcf_file, min_freq=0.0, min_depth=999999, min_qual=0) as reader:
            variants = reader.fetch("chr1", 0, 350).variants
            assert len(variants) == 0

    def test_qual_filter(self, vcf_file):
        """Variant 4 at pos 297 has QUAL=15, should be filtered at min_qual=20."""
        with VariantReader(vcf_file, min_freq=0.0, min_depth=0, min_qual=20) as reader:
            variants = reader.fetch("chr1", 230, 311).variants
            assert len(variants) == 0  # variant 4 filtered

    def test_qual_filter_disabled(self, vcf_file):
        """With min_qual=0, variant 4 should pass."""
        with VariantReader(vcf_file, min_freq=0.0, min_depth=0, min_qual=0) as reader:
            variants = reader.fetch("chr1", 230, 311).variants
            assert len(variants) == 1
            assert variants[0].pos == 296


class TestVariantAoRo:
    """Tests for ao/ro fields on Variant objects."""

    def test_pool_variant_has_ao_ro(self, vcf_file):
        """Pool-mode variants from AD field populate ao and ro."""
        with VariantReader(vcf_file, min_freq=0.0, min_depth=0, min_qual=0) as reader:
            result = reader.fetch("chr1", 0, 10)
        # pos 6 (1-based) = pos 5 (0-based): T>C, AD=80,20
        var = [v for v in result.variants if v.pos == 5][0]
        assert var.ao == 20
        assert var.ro == 80

    def test_pool_variant_ao_ro_second_variant(self, vcf_file):
        """Second variant also has correct ao/ro."""
        with VariantReader(vcf_file, min_freq=0.0, min_depth=0, min_qual=0) as reader:
            result = reader.fetch("chr1", 0, 10)
        # pos 7 (1-based) = pos 6 (0-based): G>A, AD=70,30
        var = [v for v in result.variants if v.pos == 6][0]
        assert var.ao == 30
        assert var.ro == 70

    def test_individual_variant_has_ao_ro(self, individual_vcf_file):
        """Individual-mode variants populate ao/ro from GT allele counts."""
        with IndividualVariantReader(
            individual_vcf_file, samples=["S1", "S2", "S3", "S4"],
            min_freq=0.0, min_qual=0.0, pass_only=False,
            keep_multiallelic=False, min_call_rate=0.0, min_an=0,
        ) as reader:
            result = reader.fetch("chr1", 0, 10)
        # pos 6 (1-based) = pos 5: T>C, S1:0/1 S2:0/0 S3:0/1 S4:./.
        # called=3, ref_count=4, alt_count=2
        var = [v for v in result.variants if v.pos == 5][0]
        assert var.ao == 2
        assert var.ro == 4


class TestMultiallelicFiltering:
    def test_default_skips_multiallelic(self, multiallelic_vcf_file):
        """Default behavior: positions with >1 ALT allele are discarded."""
        with VariantReader(multiallelic_vcf_file, min_freq=0.0, min_depth=0,
                           min_qual=0) as reader:
            variants = reader.fetch("chr1", 0, 350).variants
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
            variants = reader.fetch("chr1", 0, 350).variants
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
            variants = reader.fetch("chr1", 0, 350).variants
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
            variants = reader.fetch("chr1", 0, 350).variants
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
            variants = reader.fetch("chr1", 0, 350).variants
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
            variants = reader.fetch("chr1", 0, 350).variants
            assert len(variants) == 2
            positions = [v.pos for v in variants]
            assert 5 in positions
            assert 194 in positions
            assert 6 not in positions


class TestADMissingFallback:
    """Tests for _extract_freq_depth fallback paths when FORMAT/AD is absent."""

    def test_info_af_with_format_dp(self, vcf_no_ad_info_af):
        """Without AD, freq/depth come from INFO/AF + FORMAT/DP."""
        with VariantReader(vcf_no_ad_info_af, min_freq=0.0, min_depth=0,
                           min_qual=0) as reader:
            variants = reader.fetch("chr1", 0, 350).variants
            assert len(variants) == 2
            v1 = next(v for v in variants if v.pos == 5)
            assert abs(v1.freq - 0.20) < 1e-6
            assert v1.depth == 100
            v2 = next(v for v in variants if v.pos == 194)
            assert abs(v2.freq - 0.40) < 1e-6
            assert v2.depth == 100

    def test_info_af_with_info_dp_only(self, vcf_no_ad_no_format_dp):
        """Without AD and without FORMAT/DP, depth falls back to INFO/DP."""
        with VariantReader(vcf_no_ad_no_format_dp, min_freq=0.0, min_depth=0,
                           min_qual=0) as reader:
            variants = reader.fetch("chr1", 0, 350).variants
            assert len(variants) == 2
            v1 = next(v for v in variants if v.pos == 5)
            assert abs(v1.freq - 0.25) < 1e-6
            assert v1.depth == 80
            v2 = next(v for v in variants if v.pos == 194)
            assert abs(v2.freq - 0.50) < 1e-6
            assert v2.depth == 120

    def test_no_ad_no_af_returns_zero_freq(self, vcf_no_ad_no_af):
        """Without AD and without INFO/AF, freq is 0.0 — variant filtered by min_freq."""
        with VariantReader(vcf_no_ad_no_af, min_freq=0.01, min_depth=0,
                           min_qual=0) as reader:
            variants = reader.fetch("chr1", 0, 350).variants
            assert len(variants) == 0  # freq=0.0 < min_freq=0.01

    def test_no_ad_no_af_passes_with_zero_min_freq(self, vcf_no_ad_no_af):
        """Without AD/AF, freq=0.0 passes when min_freq=0."""
        with VariantReader(vcf_no_ad_no_af, min_freq=0.0, min_depth=0,
                           min_qual=0) as reader:
            variants = reader.fetch("chr1", 0, 350).variants
            assert len(variants) == 1
            assert variants[0].freq == 0.0
            assert variants[0].depth == 100

    def test_multiallelic_keep_no_ad_uses_original_freq(self, vcf_multiallelic_no_ad):
        """Multiallelic keep-mode without AD falls back to original freq/depth."""
        with VariantReader(vcf_multiallelic_no_ad, min_freq=0.0, min_depth=0,
                           min_qual=0, keep_multiallelic=True) as reader:
            variants = reader.fetch("chr1", 0, 350).variants
            # pos 7 multiallelic: two decomposed records, no AD → fallback
            ma = [v for v in variants if v.pos == 6]
            assert len(ma) == 2
            va = next(v for v in ma if v.alt == "A")
            vc = next(v for v in ma if v.alt == "C")
            assert abs(va.freq - 0.30) < 1e-6  # original INFO/AF preserved
            assert abs(vc.freq - 0.20) < 1e-6
            assert va.depth == 100
            assert vc.depth == 100

    def test_multiallelic_keep_no_ad_respects_min_freq(self, vcf_multiallelic_no_ad):
        """Multiallelic keep fallback still applies min_freq filter."""
        with VariantReader(vcf_multiallelic_no_ad, min_freq=0.25, min_depth=0,
                           min_qual=0, keep_multiallelic=True) as reader:
            variants = reader.fetch("chr1", 0, 350).variants
            ma = [v for v in variants if v.pos == 6]
            assert len(ma) == 1  # only A (0.30) passes, C (0.20) filtered
            assert ma[0].alt == "A"

    def test_multiallelic_keep_no_ad_respects_min_depth(self, vcf_multiallelic_no_ad):
        """Multiallelic keep fallback still applies min_depth filter."""
        with VariantReader(vcf_multiallelic_no_ad, min_freq=0.0, min_depth=200,
                           min_qual=0, keep_multiallelic=True) as reader:
            variants = reader.fetch("chr1", 0, 350).variants
            ma = [v for v in variants if v.pos == 6]
            assert len(ma) == 0  # depth=100 < min_depth=200

    def test_multiallelic_skip_no_ad(self, vcf_multiallelic_no_ad):
        """Default keep_multiallelic=False still skips multiallelic without AD."""
        with VariantReader(vcf_multiallelic_no_ad, min_freq=0.0, min_depth=0,
                           min_qual=0) as reader:
            variants = reader.fetch("chr1", 0, 350).variants
            positions = [v.pos for v in variants]
            assert 6 not in positions  # multiallelic skipped
            assert 5 in positions
            assert 194 in positions


class TestIndividualVariantReader:
    def test_basic_gt_frequency(self, individual_vcf_file):
        """All 4 samples, no filtering. Verify GT-derived frequencies."""
        with IndividualVariantReader(
            individual_vcf_file,
            samples=["S1", "S2", "S3", "S4"],
            min_freq=0.0, min_qual=0.0, pass_only=False,
            keep_multiallelic=False, min_call_rate=0.0, min_an=0,
        ) as reader:
            variants = reader.fetch("chr1", 0, 350).variants
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
            variants = reader.fetch("chr1", 0, 350).variants
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
            variants = reader.fetch("chr1", 0, 350).variants
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
            variants = reader.fetch("chr1", 0, 350).variants
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
            variants = reader.fetch("chr1", 0, 350).variants
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
            variants = reader.fetch("chr1", 0, 100).variants
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
            variants = reader.fetch("chr1", 0, 10).variants
            # pos 6: S4 is ./., so 0 called -> skipped
            positions = [v.pos for v in variants]
            assert 5 not in positions


    def test_duplicate_samples_uses_loaded_count(self, individual_vcf_file):
        """Duplicate sample names: n_samples should reflect loaded (deduplicated) count."""
        with IndividualVariantReader(
            individual_vcf_file, samples=["S1", "S1", "S2"],
            min_freq=0.0, min_qual=0.0, pass_only=False,
            keep_multiallelic=False, min_call_rate=0.0, min_an=0,
        ) as reader:
            # cyvcf2 deduplicates: only S1, S2 loaded
            assert reader.n_samples == 2
            variants = reader.fetch("chr1", 0, 10).variants
            v6 = next(v for v in variants if v.pos == 5)
            # S1: 0/1, S2: 0/0 -> called=2, call_rate=2/2=1.0
            assert abs(v6.call_rate - 1.0) < 1e-6

    def test_samples_none_uses_all(self, individual_vcf_file):
        """samples=None loads all VCF samples."""
        with IndividualVariantReader(
            individual_vcf_file, samples=None,
            min_freq=0.0, min_qual=0.0, pass_only=False,
            keep_multiallelic=False, min_call_rate=0.0, min_an=0,
        ) as reader:
            assert reader.n_samples == 4


class TestIndividualMultiallelic:
    def test_default_skips_multiallelic(self, individual_multiallelic_vcf_file):
        """pos 7 has ALT=A,C -> multiallelic, skipped by default."""
        with IndividualVariantReader(
            individual_multiallelic_vcf_file,
            samples=["S1", "S2", "S3", "S4"],
            min_freq=0.0, min_qual=0.0, pass_only=False,
            keep_multiallelic=False, min_call_rate=0.0, min_an=0,
        ) as reader:
            variants = reader.fetch("chr1", 0, 350).variants
            positions = [v.pos for v in variants]
            assert 6 not in positions  # multiallelic skipped
            assert 5 in positions
            assert 194 in positions
            assert len(variants) == 2

    def test_keep_multiallelic(self, individual_multiallelic_vcf_file):
        """With keep_multiallelic=True, both ALTs at pos 7 are kept.

        pos 7 G>A,C: S1:0/1 S2:0/2 S3:1/2 S4:1/2
        AN=8, AC_A=3, AC_C=3, freq_A=3/8=0.375, freq_C=3/8=0.375
        """
        with IndividualVariantReader(
            individual_multiallelic_vcf_file,
            samples=["S1", "S2", "S3", "S4"],
            min_freq=0.0, min_qual=0.0, pass_only=False,
            keep_multiallelic=True, min_call_rate=0.0, min_an=0,
        ) as reader:
            variants = reader.fetch("chr1", 0, 350).variants
            ma = [v for v in variants if v.pos == 6]
            assert len(ma) == 2
            va = next(v for v in ma if v.alt == "A")
            vc = next(v for v in ma if v.alt == "C")
            assert abs(va.freq - 3 / 8) < 1e-6
            assert abs(vc.freq - 3 / 8) < 1e-6
            assert va.depth == 8  # AN
            assert vc.depth == 8

    def test_multiallelic_min_freq_filter(self, individual_multiallelic_vcf_file):
        """min_freq filters individual ALTs at multiallelic site."""
        with IndividualVariantReader(
            individual_multiallelic_vcf_file,
            samples=["S1", "S2", "S3", "S4"],
            min_freq=0.40, min_qual=0.0, pass_only=False,
            keep_multiallelic=True, min_call_rate=0.0, min_an=0,
        ) as reader:
            variants = reader.fetch("chr1", 0, 350).variants
            # Both ALTs at pos 7 have freq=0.375 < 0.40
            ma = [v for v in variants if v.pos == 6]
            assert len(ma) == 0


class TestNonDiploidGenotypes:
    """Haploid calls are uncalled, not fabricated — see IndividualVariantReader."""

    def _vcf(self, tmp_path, name, samples, row):
        return _write_vcf(tmp_path, name, samples, [row])

    def test_all_haploid_site_is_uncalled(self, tmp_path):
        """3 hemizygous males, GT 0/0/1: no fabricated 0.167 (true AF is 1/3)."""
        path = self._vcf(tmp_path, "hap.vcf", ["m1", "m2", "m3"],
                         "chr1\t100\t.\tA\tT\t99\tPASS\t.\tGT\t0\t0\t1")
        with IndividualVariantReader(path, min_freq=0.0, min_qual=0.0,
                                     min_call_rate=0.0, min_an=1) as reader:
            assert reader.fetch("chr1", 0, 1000).variants == []

    def test_mixed_ploidy_counts_only_diploid_samples(self, tmp_path):
        """2 haploid + 1 diploid 0/1: AF from the diploid sample alone."""
        path = self._vcf(tmp_path, "mixed.vcf", ["m1", "m2", "d1"],
                         "chr1\t100\t.\tA\tT\t99\tPASS\t.\tGT\t0\t1\t0/1")
        with IndividualVariantReader(path, min_freq=0.0, min_qual=0.0,
                                     min_call_rate=0.0, min_an=1) as reader:
            (var,) = reader.fetch("chr1", 0, 1000).variants
            assert var.freq == 0.5
            assert var.depth == 2          # AN over diploid samples only
            assert var.call_rate == 1 / 3  # haploid calls lower the call rate

    def test_haploid_calls_are_counted_in_filter_stats(self, tmp_path):
        """The count rides home on FetchResult so a pool parent can warn once."""
        path = self._vcf(tmp_path, "warn.vcf", ["m1", "m2", "d1"],
                         "chr1\t100\t.\tA\tT\t99\tPASS\t.\tGT\t0\t1\t0/1")
        with IndividualVariantReader(path, min_freq=0.0, min_qual=0.0,
                                     min_call_rate=0.0, min_an=1) as reader:
            assert reader.fetch("chr1", 0, 1000).stats.n_non_diploid == 2

    def test_diploid_unaffected(self, tmp_path):
        path = self._vcf(tmp_path, "dip.vcf", ["d1", "d2", "d3"],
                         "chr1\t100\t.\tA\tT\t99\tPASS\t.\tGT\t0/0\t0/0\t0/1")
        with IndividualVariantReader(path, min_freq=0.0, min_qual=0.0,
                                     min_call_rate=0.0, min_an=1) as reader:
            (var,) = reader.fetch("chr1", 0, 1000).variants
            assert var.freq == pytest.approx(1 / 6)
            assert var.depth == 6


class TestRepeatedAlleleAtOnePosition:
    """Repeated alleles are summed — see VariantReader.fetch."""

    @staticmethod
    def _ad_vcf(tmp_path, name, rows):
        return _write_vcf(tmp_path, name, ["S1"], rows, fmt=_AD_FORMAT)

    @pytest.fixture
    def repeated_allele_vcf(self, tmp_path):
        return self._ad_vcf(tmp_path, "dup_allele.vcf", [
            "chr1\t100\t.\tG\tA\t99\tPASS\t.\tAD\t9,56",
            "chr1\t100\t.\tG\tA\t99\tPASS\t.\tAD\t9,22",
        ])

    def test_observations_are_summed(self, repeated_allele_vcf):
        with VariantReader(repeated_allele_vcf, min_freq=0.0, min_depth=0,
                           min_qual=0.0, keep_multiallelic=True) as reader:
            variants = reader.fetch("chr1", 0, 1000).variants
        assert len(variants) == 1, "one Variant per allele, not per record"
        (var,) = variants
        assert var.ao == 78    # 56 + 22
        assert var.ro == 9
        assert var.depth == 87  # ro + total alt, unchanged
        assert var.freq == pytest.approx(78 / 87)

    def test_distinct_alts_still_kept_apart(self, tmp_path):
        """Genuine multiallelic records must not be collapsed."""
        path = self._ad_vcf(tmp_path, "multi.vcf", [
            "chr1\t100\t.\tG\tA\t99\tPASS\t.\tAD\t9,56",
            "chr1\t100\t.\tG\tC\t99\tPASS\t.\tAD\t9,22",
        ])
        with VariantReader(path, min_freq=0.0, min_depth=0,
                           min_qual=0.0, keep_multiallelic=True) as reader:
            variants = reader.fetch("chr1", 0, 1000).variants
        assert sorted(v.alt for v in variants) == ["A", "C"]
        assert {v.depth for v in variants} == {87}
