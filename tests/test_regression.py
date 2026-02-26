"""Regression tests using real Acyrthosiphon pisum pool-seq dataset (400 genes, ~92k SNPs).

These tests pin numerical results validated against SNPGenie to prevent regressions
in site counting, multiallelic handling, and genome-wide summary statistics.

Default behavior skips multiallelic sites; --keep-multiallelic merges them (old default).

Data must be extracted before running:
    tar xzf data/Acyrthosiphon_pisum.tar.gz -C data/
"""

import pytest
import pandas as pd
from pathlib import Path
from click.testing import CliRunner
from pie.cli import main

REAL_DATA_DIR = Path(__file__).parents[1] / "data" / "Acyrthosiphon_pisum"


def _skip_if_no_data():
    ref = REAL_DATA_DIR / "Acyrthosiphon_pisum.fa"
    if not ref.exists():
        pytest.skip("Real test data not extracted (run: tar xzf data/Acyrthosiphon_pisum.tar.gz -C data/)")


@pytest.fixture(scope="module")
def real_results(tmp_path_factory):
    """Run pie once on the real dataset (default: skip multiallelic) and cache results."""
    _skip_if_no_data()
    ref = REAL_DATA_DIR / "Acyrthosiphon_pisum.fa"
    gff = REAL_DATA_DIR / "Acyrthosiphon_pisum.gff"
    vcf = REAL_DATA_DIR / "SRR27175631.filtered.snps.vcf.gz"

    outdir = tmp_path_factory.mktemp("real_results")
    runner = CliRunner()
    result = runner.invoke(main, [
        "run",
        "--vcf", str(vcf),
        "--gff", str(gff),
        "--fasta", str(ref),
        "--outdir", str(outdir),
        "--window-size", "0",
    ])
    assert result.exit_code == 0, f"pie run failed:\n{result.output}"

    gene_df = pd.read_csv(outdir / "gene_results.tsv", sep="\t")
    summary_df = pd.read_csv(outdir / "summary.tsv", sep="\t")
    return gene_df, summary_df


@pytest.fixture(scope="module")
def real_results_keep_multiallelic(tmp_path_factory):
    """Run pie with --keep-multiallelic and cache results."""
    _skip_if_no_data()
    ref = REAL_DATA_DIR / "Acyrthosiphon_pisum.fa"
    gff = REAL_DATA_DIR / "Acyrthosiphon_pisum.gff"
    vcf = REAL_DATA_DIR / "SRR27175631.filtered.snps.vcf.gz"

    outdir = tmp_path_factory.mktemp("real_results_keep_multi")
    runner = CliRunner()
    result = runner.invoke(main, [
        "run",
        "--vcf", str(vcf),
        "--gff", str(gff),
        "--fasta", str(ref),
        "--outdir", str(outdir),
        "--window-size", "0",
        "--keep-multiallelic",
    ])
    assert result.exit_code == 0, f"pie run --keep-multiallelic failed:\n{result.output}"

    gene_df = pd.read_csv(outdir / "gene_results.tsv", sep="\t")
    summary_df = pd.read_csv(outdir / "summary.tsv", sep="\t")
    return gene_df, summary_df


@pytest.mark.slow
class TestRealDataRegression:

    def test_all_400_genes_processed(self, real_results):
        """All 400 genes across 4 chromosomes (100 each) are in the output."""
        gene_df, _ = real_results
        assert len(gene_df) == 400

        per_chrom = gene_df["chrom"].value_counts()
        for chrom in ["1", "2", "3", "X"]:
            assert per_chrom[chrom] == 100, f"Expected 100 genes on chrom {chrom}"

    def test_genome_wide_summary(self, real_results):
        """Pin genome-wide piN, piS, and piN/piS (default: skip multiallelic, exclude stops)."""
        _, summary_df = real_results
        row = summary_df.iloc[0]

        assert row["total_genes"] == 400
        assert abs(row["genome_piN"] - 0.000974) < 1e-6
        assert abs(row["genome_piS"] - 0.003939) < 1e-6
        assert abs(row["genome_piN_piS"] - 0.2472) < 1e-4

    def test_sites_sum_leq_three_per_codon(self, real_results):
        """(N_sites + S_sites) / n_codons <= 3.0 for every gene.

        With stop_gained excluded (default), mutations to stop codons reduce
        the effective number of sites at affected positions, so the per-codon
        average may be slightly below 3.0.
        """
        gene_df, _ = real_results
        sites_per_codon = (gene_df["N_sites"] + gene_df["S_sites"]) / gene_df["n_codons"]
        assert (sites_per_codon <= 3.0 + 1e-10).all()
        assert (sites_per_codon > 2.9).all()

    def test_no_negative_values(self, real_results):
        """All site counts, diffs, and diversity values must be non-negative."""
        gene_df, _ = real_results
        for col in ["N_sites", "S_sites", "N_diffs", "S_diffs", "piN", "piS"]:
            assert (gene_df[col] >= 0).all(), f"Negative values in {col}"

    def test_multiallelic_gene_apisum_017038(self, real_results):
        """Pin Apisum_017038 (default: multiallelic sites skipped, + strand)."""
        gene_df, _ = real_results
        g = gene_df[gene_df["gene_id"].str.contains("Apisum_017038")].iloc[0]

        assert g["n_codons"] == 619
        assert g["n_variants"] == 41
        assert abs(g["N_sites"] - 1421.654100) < 1e-4
        assert abs(g["S_sites"] - 435.345900) < 1e-4
        assert abs(g["piN"] - 0.006868) < 1e-5
        assert abs(g["piS"] - 0.014170) < 1e-5
        assert abs(g["piN_piS"] - 0.484712) < 1e-4

    def test_multiallelic_gene_apisum_003665(self, real_results):
        """Pin Apisum_003665 (default: multiallelic sites skipped, - strand)."""
        gene_df, _ = real_results
        g = gene_df[gene_df["gene_id"].str.contains("Apisum_003665")].iloc[0]

        assert g["n_codons"] == 215
        assert g["n_variants"] == 18
        assert abs(g["N_sites"] - 495.546269) < 1e-4
        assert abs(g["S_sites"] - 149.453731) < 1e-4
        assert abs(g["piN"] - 0.007547) < 1e-5
        assert abs(g["piS"] - 0.013266) < 1e-5
        assert abs(g["piN_piS"] - 0.568877) < 1e-4

    def test_high_variant_gene(self, real_results):
        """Pin Apisum_003662 (default: multiallelic sites skipped)."""
        gene_df, _ = real_results
        g = gene_df[gene_df["gene_id"].str.contains("Apisum_003662")].iloc[0]

        assert g["n_codons"] == 418
        assert g["n_variants"] == 79
        assert g["n_poly_codons"] == 72
        assert abs(g["piN"] - 0.012722) < 1e-5
        assert abs(g["piS"] - 0.013534) < 1e-5
        assert abs(g["piN_piS"] - 0.939993) < 1e-4

    def test_both_strands_present(self, real_results):
        """Both + and - strand genes have non-zero piN values."""
        gene_df, _ = real_results
        plus_genes = gene_df[gene_df["strand"] == "+"]
        minus_genes = gene_df[gene_df["strand"] == "-"]

        assert len(plus_genes) > 0
        assert len(minus_genes) > 0
        assert (plus_genes["piN"] > 0).any(), "No + strand genes with piN > 0"
        assert (minus_genes["piN"] > 0).any(), "No - strand genes with piN > 0"


@pytest.mark.slow
class TestRealDataKeepMultiallelic:
    """Verify --keep-multiallelic reproduces the old (merge) behavior."""

    def test_genome_wide_summary_keep_multiallelic(self, real_results_keep_multiallelic):
        """Pin genome-wide values with --keep-multiallelic."""
        _, summary_df = real_results_keep_multiallelic
        row = summary_df.iloc[0]

        assert row["total_genes"] == 400
        assert abs(row["genome_piN"] - 0.001047) < 1e-6
        assert abs(row["genome_piS"] - 0.004088) < 1e-6
        assert abs(row["genome_piN_piS"] - 0.2561) < 1e-4

    def test_multiallelic_gene_apisum_017038_keep(self, real_results_keep_multiallelic):
        """Pin Apisum_017038 with --keep-multiallelic (+ strand)."""
        gene_df, _ = real_results_keep_multiallelic
        g = gene_df[gene_df["gene_id"].str.contains("Apisum_017038")].iloc[0]

        assert g["n_codons"] == 619
        assert g["n_variants"] == 45
        assert abs(g["N_sites"] - 1421.634246) < 1e-4
        assert abs(g["S_sites"] - 435.365754) < 1e-4
        assert abs(g["piN"] - 0.007146) < 1e-5
        assert abs(g["piS"] - 0.014169) < 1e-5
        assert abs(g["piN_piS"] - 0.504338) < 1e-4

    def test_multiallelic_gene_apisum_003665_keep(self, real_results_keep_multiallelic):
        """Pin Apisum_003665 with --keep-multiallelic (- strand)."""
        gene_df, _ = real_results_keep_multiallelic
        g = gene_df[gene_df["gene_id"].str.contains("Apisum_003665")].iloc[0]

        assert g["n_codons"] == 215
        assert g["n_variants"] == 22
        assert abs(g["N_sites"] - 495.546269) < 1e-4
        assert abs(g["S_sites"] - 149.453731) < 1e-4
        assert abs(g["piN"] - 0.008552) < 1e-5
        assert abs(g["piS"] - 0.017477) < 1e-5
        assert abs(g["piN_piS"] - 0.489328) < 1e-4

    def test_high_variant_gene_keep(self, real_results_keep_multiallelic):
        """Pin Apisum_003662 with --keep-multiallelic."""
        gene_df, _ = real_results_keep_multiallelic
        g = gene_df[gene_df["gene_id"].str.contains("Apisum_003662")].iloc[0]

        assert g["n_codons"] == 418
        assert g["n_variants"] == 93
        assert g["n_poly_codons"] == 77
        assert abs(g["piN"] - 0.014266) < 1e-5
        assert abs(g["piS"] - 0.013886) < 1e-5
        assert abs(g["piN_piS"] - 1.027330) < 1e-4

    def test_more_variants_with_keep_multiallelic(self, real_results, real_results_keep_multiallelic):
        """--keep-multiallelic should have >= variants than default for affected genes."""
        gene_df_default, _ = real_results
        gene_df_keep, _ = real_results_keep_multiallelic

        total_default = gene_df_default["n_variants"].sum()
        total_keep = gene_df_keep["n_variants"].sum()
        assert total_keep > total_default, \
            f"Expected more variants with --keep-multiallelic ({total_keep}) than default ({total_default})"
