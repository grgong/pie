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
        "pool",
        "--vcf", str(vcf),
        "--gff", str(gff),
        "--fasta", str(ref),
        "--outdir", str(outdir),
        "--window-size", "0",
    ])
    assert result.exit_code == 0, f"pie pool failed:\n{result.output}"

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
        "pool",
        "--vcf", str(vcf),
        "--gff", str(gff),
        "--fasta", str(ref),
        "--outdir", str(outdir),
        "--window-size", "0",
        "--keep-multiallelic",
    ])
    assert result.exit_code == 0, f"pie pool --keep-multiallelic failed:\n{result.output}"

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
        assert abs(g["N_sites"] - 1420.987433) < 1e-4
        assert abs(g["S_sites"] - 436.012567) < 1e-4
        assert abs(g["piN"] - 0.006871) < 1e-5
        assert abs(g["piS"] - 0.014148) < 1e-5
        assert abs(g["piN_piS"] - 0.485682) < 1e-4

    def test_multiallelic_gene_apisum_003665(self, real_results):
        """Pin Apisum_003665 (default: multiallelic sites skipped, - strand)."""
        gene_df, _ = real_results
        g = gene_df[gene_df["gene_id"].str.contains("Apisum_003665")].iloc[0]

        assert g["n_codons"] == 215
        assert g["n_variants"] == 18
        assert abs(g["N_sites"] - 494.879603) < 1e-4
        assert abs(g["S_sites"] - 150.120397) < 1e-4
        assert abs(g["piN"] - 0.007557) < 1e-5
        assert abs(g["piS"] - 0.013207) < 1e-5
        assert abs(g["piN_piS"] - 0.572184) < 1e-4

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
    """Pin --keep-multiallelic output.

    Values were re-pinned when allele observations stopped being overwritten
    across records naming the same allele (haplotype-decomposed VCF): see
    TestHaplotypeDecomposedRecords below. Only --keep-multiallelic runs are
    affected — the default mode drops these sites as multiallelic.
    """

    def test_genome_wide_summary_keep_multiallelic(self, real_results_keep_multiallelic):
        """Pin genome-wide values with --keep-multiallelic."""
        _, summary_df = real_results_keep_multiallelic
        row = summary_df.iloc[0]

        assert row["total_genes"] == 400
        assert abs(row["genome_piN"] - 0.001051) < 1e-6
        assert abs(row["genome_piS"] - 0.004082) < 1e-6
        assert abs(row["genome_piN_piS"] - 0.2574) < 1e-4

    def test_multiallelic_gene_apisum_017038_keep(self, real_results_keep_multiallelic):
        """Pin Apisum_017038 with --keep-multiallelic (+ strand)."""
        gene_df, _ = real_results_keep_multiallelic
        g = gene_df[gene_df["gene_id"].str.contains("Apisum_017038")].iloc[0]

        assert g["n_codons"] == 619
        assert g["n_variants"] == 43
        assert abs(g["N_sites"] - 1420.291775) < 1e-4
        assert abs(g["S_sites"] - 436.708225) < 1e-4
        assert abs(g["piN"] - 0.007211) < 1e-5
        assert abs(g["piS"] - 0.014125) < 1e-5
        assert abs(g["piN_piS"] - 0.510521) < 1e-4

    def test_multiallelic_gene_apisum_003665_keep(self, real_results_keep_multiallelic):
        """Pin Apisum_003665 with --keep-multiallelic (- strand)."""
        gene_df, _ = real_results_keep_multiallelic
        g = gene_df[gene_df["gene_id"].str.contains("Apisum_003665")].iloc[0]

        assert g["n_codons"] == 215
        assert g["n_variants"] == 21
        assert abs(g["N_sites"] - 494.879603) < 1e-4
        assert abs(g["S_sites"] - 150.120397) < 1e-4
        assert abs(g["piN"] - 0.008274) < 1e-5
        assert abs(g["piS"] - 0.017400) < 1e-5
        assert abs(g["piN_piS"] - 0.475538) < 1e-4

    def test_high_variant_gene_keep(self, real_results_keep_multiallelic):
        """Pin Apisum_003662 with --keep-multiallelic."""
        gene_df, _ = real_results_keep_multiallelic
        g = gene_df[gene_df["gene_id"].str.contains("Apisum_003662")].iloc[0]

        assert g["n_codons"] == 418
        assert g["n_variants"] == 87
        assert g["n_poly_codons"] == 77
        assert abs(g["piN"] - 0.013910) < 1e-5
        assert abs(g["piS"] - 0.013879) < 1e-5
        assert abs(g["piN_piS"] - 1.002267) < 1e-4

    def test_more_variants_with_keep_multiallelic(self, real_results, real_results_keep_multiallelic):
        """--keep-multiallelic should have >= variants than default for affected genes."""
        gene_df_default, _ = real_results
        gene_df_keep, _ = real_results_keep_multiallelic

        total_default = gene_df_default["n_variants"].sum()
        total_keep = gene_df_keep["n_variants"].sum()
        assert total_keep > total_default, \
            f"Expected more variants with --keep-multiallelic ({total_keep}) than default ({total_default})"


@pytest.fixture(scope="module")
def real_variant_table(tmp_path_factory):
    """Run pie once with --keep-multiallelic --variant-table; cache the table."""
    _skip_if_no_data()
    ref = REAL_DATA_DIR / "Acyrthosiphon_pisum.fa"
    gff = REAL_DATA_DIR / "Acyrthosiphon_pisum.gff"
    vcf = REAL_DATA_DIR / "SRR27175631.filtered.snps.vcf.gz"

    outdir = tmp_path_factory.mktemp("real_variant_table")
    result = CliRunner().invoke(main, [
        "pool", "--vcf", str(vcf), "--gff", str(gff), "--fasta", str(ref),
        "--outdir", str(outdir), "--keep-multiallelic", "--variant-table",
    ])
    assert result.exit_code == 0, result.output
    return pd.read_csv(outdir / "variant_results.tsv", sep="\t")


@pytest.mark.slow
class TestHaplotypeDecomposedRecords:
    """Repeated alleles are summed — see VariantReader.fetch.

    chr3:1714 in this dataset is `G>A` twice, AO 56 and 22, RO 9: one allele
    on two haplotype backgrounds, so the site's A frequency is 78/87.
    """

    def test_repeated_allele_is_summed(self, real_variant_table):
        site = real_variant_table[
            (real_variant_table["chrom"].astype(str) == "3")
            & (real_variant_table["pos"] == 1714)]
        assert len(site) == 1, "one row per allele, not one per VCF record"
        row = site.iloc[0]
        assert row["ao"] == 78  # 56 + 22
        assert row["ro"] == 9
        assert row["dp"] == 87
        assert abs(row["af"] - 78 / 87) < 1e-6

    def test_no_allele_appears_twice_at_a_position(self, real_variant_table):
        counts = real_variant_table.groupby(
            ["chrom", "pos", "ref", "alt"]).size()
        assert (counts == 1).all(), f"{(counts > 1).sum()} duplicated alleles"
