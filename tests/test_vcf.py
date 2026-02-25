import pytest
import os
import shutil
from pie.vcf import ensure_indexed, VariantReader


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
