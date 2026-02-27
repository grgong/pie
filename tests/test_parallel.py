"""Tests for gene-level multiprocessing parallel runner."""

import pytest
from pie.annotation import NoGenesFoundError
from pie.parallel import run_parallel
from pie.diversity import GeneResult


class TestRunParallel:
    def test_single_thread(self, ref_fasta, gff3_file, vcf_file):
        results = run_parallel(ref_fasta, gff3_file, vcf_file,
                               min_freq=0.0, min_depth=0, min_qual=0,
                               pass_only=False, threads=1)
        assert len(results) == 3
        assert all(isinstance(r, GeneResult) for r in results)

    def test_multi_thread(self, ref_fasta, gff3_file, vcf_file):
        results = run_parallel(ref_fasta, gff3_file, vcf_file,
                               min_freq=0.0, min_depth=0, min_qual=0,
                               pass_only=False, threads=2)
        assert len(results) == 3

    def test_results_match(self, ref_fasta, gff3_file, vcf_file):
        r1 = run_parallel(ref_fasta, gff3_file, vcf_file,
                          min_freq=0.0, min_depth=0, min_qual=0,
                          pass_only=False, threads=1)
        r2 = run_parallel(ref_fasta, gff3_file, vcf_file,
                          min_freq=0.0, min_depth=0, min_qual=0,
                          pass_only=False, threads=2)
        for a, b in zip(r1, r2):
            assert a.gene_id == b.gene_id
            assert abs(a.piN - b.piN) < 1e-10
            assert abs(a.piS - b.piS) < 1e-10

    def test_sorted_by_position(self, ref_fasta, gff3_file, vcf_file):
        results = run_parallel(ref_fasta, gff3_file, vcf_file,
                               min_freq=0.0, min_depth=0, min_qual=0,
                               threads=1)
        positions = [(r.chrom, r.start) for r in results]
        assert positions == sorted(positions)

    def test_with_qual_filter(self, ref_fasta, gff3_file, vcf_file):
        results = run_parallel(ref_fasta, gff3_file, vcf_file,
                               min_freq=0.0, min_depth=0, min_qual=20,
                               threads=1)
        gene3 = [r for r in results if "gene3" in r.gene_id.lower()][0]
        assert gene3.n_variants == 0  # QUAL=15 variant filtered


class TestNoGenesFoundError:
    """Issue #4: CDS-only GFF yields 0 genes -> NoGenesFoundError."""

    def test_cdsonly_gff_raises(self, ref_fasta, cdsonly_gff, vcf_file):
        with pytest.raises(NoGenesFoundError, match="No genes with CDS features"):
            run_parallel(ref_fasta, cdsonly_gff, vcf_file,
                         min_freq=0.0, min_depth=0, min_qual=0, threads=1)


class TestContigMismatch:
    """Issue #7: Contig naming mismatch detection."""

    def test_chr1_vs_1_raises(self, ref_fasta, gff3_file, mismatch_vcf_file):
        """GFF has 'chr1', VCF has '1' -> ValueError with helpful message."""
        with pytest.raises(ValueError, match="No contig names shared"):
            run_parallel(ref_fasta, gff3_file, mismatch_vcf_file,
                         min_freq=0.0, min_depth=0, min_qual=0, threads=1)
