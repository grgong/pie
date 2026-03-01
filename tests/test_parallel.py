"""Tests for gene-level multiprocessing parallel runner."""

import logging
import pytest
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

    def test_stop_codon_summary_logged_once(self, ref_fasta, gff3_file, vcf_file, caplog):
        """Stop-codon renormalization is summarized once, not per-gene."""
        with caplog.at_level(logging.WARNING, logger="pie.parallel"):
            results = run_parallel(ref_fasta, gff3_file, vcf_file,
                                   min_freq=0.0, min_depth=0, min_qual=0,
                                   threads=1)
        stop_msgs = [r for r in caplog.records
                     if "Stop-codon renormalization" in r.message]
        # Either 0 (no genes had stop codons) or exactly 1 summary line
        assert len(stop_msgs) <= 1

    def test_n_stop_codons_on_results(self, ref_fasta, gff3_file, vcf_file):
        """GeneResult carries n_stop_codons field."""
        results = run_parallel(ref_fasta, gff3_file, vcf_file,
                               min_freq=0.0, min_depth=0, min_qual=0,
                               threads=1)
        for r in results:
            assert hasattr(r, "n_stop_codons")
            assert r.n_stop_codons >= 0
