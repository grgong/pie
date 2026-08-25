"""Tests for gene-level multiprocessing parallel runner."""

import logging
import pytest
from pie.parallel import run_parallel
from pie.diversity import GeneResult
from tests.helpers import bgzip_and_index


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


class TestWorkerInitFailure:
    """Bad inputs must raise, not hang — see run_parallel's pre-flight."""

    def test_bad_fasta_raises_instead_of_hanging(self, gff3_file, vcf_file):
        for threads in (1, 2):
            with pytest.raises(OSError):
                run_parallel("/nonexistent/ref.fa", gff3_file, vcf_file,
                             threads=threads)


class TestAtexitRegistration:
    """No atexit hook anywhere — a pool worker would never run one anyway."""

    def test_run_registers_no_atexit_hook(self, ref_fasta, gff3_file,
                                          vcf_file):
        import atexit

        before = atexit._ncallbacks()
        for _ in range(3):
            run_parallel(ref_fasta, gff3_file, vcf_file,
                         min_freq=0.0, min_depth=0, min_qual=0, threads=1)
        assert atexit._ncallbacks() == before


@pytest.fixture
def haploid_vcf_file(tmp_path):
    """Individual-mode VCF whose pos-195 site carries two haploid calls."""
    vcf_content = (
        "##fileformat=VCFv4.2\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "##contig=<ID=chr1,length=350>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\n"
        "chr1\t195\t.\tA\tT\t45\t.\t.\tGT\t0/1\t0/1\t1\t0\n"
    )
    vcf_path = tmp_path / "haploid.vcf"
    vcf_path.write_text(vcf_content)
    return bgzip_and_index(vcf_path)


class TestReaderDiagnostics:
    """Reader observations must reach the user in pool runs too.

    A worker's reader is never closed by anyone — the pool terminates it — so
    these counts travel back in FilterStats and are summed by run_parallel.
    """

    @pytest.mark.parametrize("threads", [1, 2])
    def test_non_diploid_calls_warned_once(self, ref_fasta, gff3_file,
                                           haploid_vcf_file, caplog, threads):
        with caplog.at_level(logging.WARNING, logger="pie.parallel"):
            run_parallel(ref_fasta, gff3_file, haploid_vcf_file,
                         min_freq=0.0, min_qual=0.0, mode="individual",
                         min_call_rate=0.0, min_an=1, threads=threads)
        msgs = [r.message for r in caplog.records
                if "non-diploid genotype call(s)" in r.message]
        assert len(msgs) == 1
        assert "2 non-diploid" in msgs[0]

    def test_silent_when_all_diploid(self, ref_fasta, gff3_file,
                                     individual_vcf_file, caplog):
        with caplog.at_level(logging.WARNING, logger="pie.parallel"):
            run_parallel(ref_fasta, gff3_file, individual_vcf_file,
                         min_freq=0.0, min_qual=0.0, mode="individual",
                         min_call_rate=0.0, min_an=1, threads=2)
        assert "non-diploid" not in caplog.text
