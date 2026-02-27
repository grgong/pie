"""Integration tests verifying complete pipeline against hand-calculated piN/piS values.

These tests run the full CLI pipeline end-to-end and compare the output TSV
values with manually verified Nei-Gojobori calculations.  See
tests/create_test_data.py for the detailed derivation of all expected numbers.

Test data layout (chr1, 350 bp):
  Gene1: pos 1-90,   + strand, 30 codons (29 excl. stop), 2 variants
  Gene2: pos 101-220, + strand, 2 exons (boundary-spanning codon), 1 variant
  Gene3: pos 231-311, - strand, 27 codons (26 excl. stop), 1 low-QUAL variant
"""

import pytest
import pandas as pd
from click.testing import CliRunner
from pie.cli import main


@pytest.fixture
def runner():
    return CliRunner()


class TestEndToEnd:
    def _run_analysis(self, runner, ref_fasta, gff3_file, vcf_file, tmp_path, min_qual=0):
        result = runner.invoke(main, [
            "run",
            "--vcf", vcf_file,
            "--gff", gff3_file,
            "--fasta", ref_fasta,
            "--outdir", str(tmp_path),
            "--min-freq", "0",
            "--min-depth", "0",
            "--min-qual", str(min_qual),
            "--window-size", "30",
            "--window-step", "10",
        ])
        assert result.exit_code == 0, result.output
        return tmp_path

    def test_gene1_piN_piS(self, runner, ref_fasta, gff3_file, vcf_file, tmp_path):
        """Gene1: ATG GCT GAT GCT*26 TAA (29 coding codons)
        Variant 1: pos 6 T->C, GCT->GCC (syn), freq=0.20
        Variant 2: pos 7 G->A, GAT->AAT (nonsyn), freq=0.30

        N_sites = 3.0 + 8/3 + 27*2.0 = 59.6667
        S_sites = 0.0 + 1/3 + 27*1.0 = 27.3333
        N_diffs = 2*0.70*0.30 = 0.42
        S_diffs = 2*0.80*0.20 = 0.32
        piN = 0.42/59.6667 = 0.007039
        piS = 0.32/27.3333 = 0.011707
        piN/piS = 0.6013
        """
        out = self._run_analysis(runner, ref_fasta, gff3_file, vcf_file, tmp_path)
        df = pd.read_csv(out / "gene_results.tsv", sep="\t")
        g1 = df[df["gene_id"].str.contains("gene1", case=False)].iloc[0]

        assert g1["n_codons"] == 29
        assert g1["n_poly_codons"] == 2
        assert g1["n_variants"] == 2
        assert abs(g1["N_sites"] - 59.6667) < 0.001
        assert abs(g1["S_sites"] - 27.3333) < 0.001
        assert abs(g1["N_diffs"] - 0.42) < 1e-6
        assert abs(g1["S_diffs"] - 0.32) < 1e-6
        assert abs(g1["piN"] - 0.007039) < 0.0001
        assert abs(g1["piS"] - 0.011707) < 0.0001

    def test_gene2_multi_exon(self, runner, ref_fasta, gff3_file, vcf_file, tmp_path):
        """Gene2: 2 exons, 32 coding codons, 1 nonsynonymous variant
        Variant: pos 195 A->T, GAT->GTT (Asp->Val), freq=0.40

        N_diffs = 2*0.60*0.40 = 0.48
        S_diffs = 0.0
        piS = 0 -> piN_piS = NA
        """
        out = self._run_analysis(runner, ref_fasta, gff3_file, vcf_file, tmp_path)
        df = pd.read_csv(out / "gene_results.tsv", sep="\t")
        g2 = df[df["gene_id"].str.contains("gene2", case=False)].iloc[0]

        assert g2["n_codons"] == 32
        assert g2["n_poly_codons"] == 1
        assert abs(g2["N_diffs"] - 0.48) < 1e-6
        assert g2["S_diffs"] == 0.0
        assert pd.isna(g2["piN_piS"])

    def test_gene3_qual_filter(self, runner, ref_fasta, gff3_file, vcf_file, tmp_path):
        """Gene3: - strand, variant QUAL=15 < 20 -> filtered with min_qual=20."""
        out = self._run_analysis(runner, ref_fasta, gff3_file, vcf_file, tmp_path, min_qual=20)
        df = pd.read_csv(out / "gene_results.tsv", sep="\t")
        g3 = df[df["gene_id"].str.contains("gene3", case=False)].iloc[0]

        assert g3["n_variants"] == 0
        assert g3["piN"] == 0.0
        assert g3["piS"] == 0.0

    def test_gene3_no_qual_filter(self, runner, ref_fasta, gff3_file, vcf_file, tmp_path):
        """Gene3 without QUAL filter: synonymous variant included.
        Variant: pos 297 A->G (fwd), T->C in CDS sense, GAT->GAC (Asp->Asp), freq=0.50
        S_diffs = 2*0.50*0.50 = 0.50
        """
        out = self._run_analysis(runner, ref_fasta, gff3_file, vcf_file, tmp_path, min_qual=0)
        df = pd.read_csv(out / "gene_results.tsv", sep="\t")
        g3 = df[df["gene_id"].str.contains("gene3", case=False)].iloc[0]

        assert g3["n_variants"] == 1
        assert abs(g3["S_diffs"] - 0.50) < 1e-6
        assert g3["N_diffs"] == 0.0

    def test_all_output_files_created(self, runner, ref_fasta, gff3_file, vcf_file, tmp_path):
        out = self._run_analysis(runner, ref_fasta, gff3_file, vcf_file, tmp_path)
        assert (out / "gene_results.tsv").exists()
        assert (out / "window_results.tsv").exists()
        assert (out / "summary.tsv").exists()

    def test_summary_genome_wide(self, runner, ref_fasta, gff3_file, vcf_file, tmp_path):
        out = self._run_analysis(runner, ref_fasta, gff3_file, vcf_file, tmp_path)
        df = pd.read_csv(out / "summary.tsv", sep="\t")
        assert len(df) == 1
        row = df.iloc[0]
        assert row["total_genes"] == 3
        # genome_piN = (0.42+0.48+0) / (59.6667+65.6667+N_sites_gene3)
        assert row["genome_piN"] > 0

    def test_window_results_have_content(self, runner, ref_fasta, gff3_file, vcf_file, tmp_path):
        out = self._run_analysis(runner, ref_fasta, gff3_file, vcf_file, tmp_path)
        df = pd.read_csv(out / "window_results.tsv", sep="\t")
        assert len(df) > 0
        assert "piN" in df.columns
        assert "piS" in df.columns

    def test_plot_generation(self, runner, ref_fasta, gff3_file, vcf_file, tmp_path):
        out = self._run_analysis(runner, ref_fasta, gff3_file, vcf_file, tmp_path)
        result = runner.invoke(main, [
            "plot",
            "--gene-results", str(out / "gene_results.tsv"),
            "--output", str(out / "manhattan.png"),
        ])
        assert result.exit_code == 0
        assert (out / "manhattan.png").exists()

    def test_keep_multiallelic_flag(self, runner, ref_fasta, gff3_file, vcf_file, tmp_path):
        """--keep-multiallelic flag is accepted and does not break pipeline."""
        result = runner.invoke(main, [
            "run",
            "--vcf", vcf_file,
            "--gff", gff3_file,
            "--fasta", ref_fasta,
            "--outdir", str(tmp_path),
            "--min-freq", "0",
            "--min-depth", "0",
            "--min-qual", "0",
            "--keep-multiallelic",
        ])
        assert result.exit_code == 0, result.output
        assert (tmp_path / "gene_results.tsv").exists()

    def test_gtf_also_works(self, runner, ref_fasta, gtf_file, vcf_file, tmp_path):
        """GTF should produce same results as GFF3."""
        result = runner.invoke(main, [
            "run", "--vcf", vcf_file, "--gff", gtf_file, "--fasta", ref_fasta,
            "--outdir", str(tmp_path), "--min-freq", "0", "--min-depth", "0", "--min-qual", "0",
        ])
        assert result.exit_code == 0
        df = pd.read_csv(tmp_path / "gene_results.tsv", sep="\t")
        assert len(df) == 3


class TestCrossFormat:
    """GFF3 and GTF should produce identical numerical results."""

    def test_gff3_gtf_agreement(self, runner, ref_fasta, gff3_file, gtf_file, vcf_file, tmp_path):
        gff_out = tmp_path / "gff3"
        gtf_out = tmp_path / "gtf"
        gff_out.mkdir()
        gtf_out.mkdir()

        for gff, outdir in [(gff3_file, gff_out), (gtf_file, gtf_out)]:
            result = runner.invoke(main, [
                "run", "--vcf", vcf_file, "--gff", gff, "--fasta", ref_fasta,
                "--outdir", str(outdir), "--min-freq", "0", "--min-depth", "0",
                "--min-qual", "0",
            ])
            assert result.exit_code == 0, result.output

        df_gff = pd.read_csv(gff_out / "gene_results.tsv", sep="\t")
        df_gtf = pd.read_csv(gtf_out / "gene_results.tsv", sep="\t")

        assert len(df_gff) == len(df_gtf)
        for col in ["N_sites", "S_sites", "N_diffs", "S_diffs", "piN", "piS"]:
            for i in range(len(df_gff)):
                assert abs(df_gff.iloc[i][col] - df_gtf.iloc[i][col]) < 1e-6, \
                    f"Mismatch in {col} for gene {df_gff.iloc[i]['gene_id']}"


class TestIndividualMode:
    """End-to-end tests for --mode individual using GT-derived frequencies."""

    def _run_individual(self, runner, ref_fasta, gff3_file, individual_vcf_file,
                        tmp_path, min_qual=0, min_call_rate=0):
        result = runner.invoke(main, [
            "run", "--mode", "individual",
            "--vcf", individual_vcf_file,
            "--gff", gff3_file,
            "--fasta", ref_fasta,
            "--outdir", str(tmp_path),
            "--min-freq", "0",
            "--min-qual", str(min_qual),
            "--min-call-rate", str(min_call_rate),
            "--min-an", "0",
            "--window-size", "30",
            "--window-step", "10",
        ])
        assert result.exit_code == 0, result.output
        return tmp_path

    def test_gene1_individual(self, runner, ref_fasta, gff3_file,
                               individual_vcf_file, tmp_path):
        """Gene1: pos6 T>C syn freq=1/3, pos7 G>A nonsyn freq=1/4."""
        out = self._run_individual(runner, ref_fasta, gff3_file,
                                    individual_vcf_file, tmp_path)
        df = pd.read_csv(out / "gene_results.tsv", sep="\t")
        g1 = df[df["gene_id"].str.contains("gene1", case=False)].iloc[0]

        assert g1["n_codons"] == 29
        assert g1["n_poly_codons"] == 2
        assert g1["n_variants"] == 2
        assert abs(g1["N_sites"] - 59.6667) < 0.001
        assert abs(g1["S_sites"] - 27.3333) < 0.001
        assert abs(g1["N_diffs"] - 0.375) < 1e-6    # 2*(3/4)*(1/4)
        assert abs(g1["S_diffs"] - 4 / 9) < 1e-6    # 2*(2/3)*(1/3)

    def test_gene2_individual(self, runner, ref_fasta, gff3_file,
                               individual_vcf_file, tmp_path):
        """Gene2: pos195 A>T nonsyn freq=5/8."""
        out = self._run_individual(runner, ref_fasta, gff3_file,
                                    individual_vcf_file, tmp_path)
        df = pd.read_csv(out / "gene_results.tsv", sep="\t")
        g2 = df[df["gene_id"].str.contains("gene2", case=False)].iloc[0]

        assert g2["n_poly_codons"] == 1
        assert abs(g2["N_diffs"] - 15 / 32) < 1e-6  # 2*(3/8)*(5/8)
        assert g2["S_diffs"] == 0.0

    def test_gene3_individual(self, runner, ref_fasta, gff3_file,
                               individual_vcf_file, tmp_path):
        """Gene3: pos297 syn freq=1/4."""
        out = self._run_individual(runner, ref_fasta, gff3_file,
                                    individual_vcf_file, tmp_path)
        df = pd.read_csv(out / "gene_results.tsv", sep="\t")
        g3 = df[df["gene_id"].str.contains("gene3", case=False)].iloc[0]

        assert abs(g3["S_diffs"] - 0.375) < 1e-6    # 2*(3/4)*(1/4)
        assert g3["N_diffs"] == 0.0

    def test_output_has_sample_metadata(self, runner, ref_fasta, gff3_file,
                                         individual_vcf_file, tmp_path):
        """Individual mode outputs include n_samples and mean_call_rate."""
        out = self._run_individual(runner, ref_fasta, gff3_file,
                                    individual_vcf_file, tmp_path)
        df = pd.read_csv(out / "gene_results.tsv", sep="\t")
        assert "n_samples" in df.columns
        assert "mean_call_rate" in df.columns
        assert (df["n_samples"] == 4).all()

    def test_summary_has_sample_metadata(self, runner, ref_fasta, gff3_file,
                                          individual_vcf_file, tmp_path):
        out = self._run_individual(runner, ref_fasta, gff3_file,
                                    individual_vcf_file, tmp_path)
        df = pd.read_csv(out / "summary.tsv", sep="\t")
        assert "n_samples_selected" in df.columns
        assert df.iloc[0]["n_samples_selected"] == 4

    def test_min_call_rate_filters_variants(self, runner, ref_fasta, gff3_file,
                                             individual_vcf_file, tmp_path):
        """min_call_rate=0.8 filters pos6 (0.75) and pos297 (0.50)."""
        out = self._run_individual(runner, ref_fasta, gff3_file,
                                    individual_vcf_file, tmp_path,
                                    min_call_rate=0.8)
        df = pd.read_csv(out / "gene_results.tsv", sep="\t")
        g1 = df[df["gene_id"].str.contains("gene1", case=False)].iloc[0]
        # Only pos7 passes -> 1 variant in gene1
        assert g1["n_variants"] == 1

    def test_all_output_files_created(self, runner, ref_fasta, gff3_file,
                                       individual_vcf_file, tmp_path):
        out = self._run_individual(runner, ref_fasta, gff3_file,
                                    individual_vcf_file, tmp_path)
        assert (out / "gene_results.tsv").exists()
        assert (out / "window_results.tsv").exists()
        assert (out / "summary.tsv").exists()


class TestRobustness:
    """Integration tests for PR#9 robustness fixes (Issues #2, #4, #7)."""

    def test_ambiguous_bases_complete(self, runner, ref_with_n_fasta, gff3_file,
                                      vcf_file, tmp_path):
        """Issue #2: N bases in reference genome should not crash (KeyError).

        Gene1 codon 2 has N -> skipped, pipeline continues with 28 codons.
        """
        result = runner.invoke(main, [
            "run", "--vcf", vcf_file, "--gff", gff3_file,
            "--fasta", ref_with_n_fasta, "--outdir", str(tmp_path),
            "--min-freq", "0", "--min-depth", "0", "--min-qual", "0",
        ])
        assert result.exit_code == 0, result.output
        df = pd.read_csv(tmp_path / "gene_results.tsv", sep="\t")
        assert len(df) == 3
        g1 = df[df["gene_id"].str.contains("gene1", case=False)].iloc[0]
        assert g1["n_codons"] == 28  # 29 - 1 skipped N codon

    def test_all_n_gene_completes(self, runner, ref_all_n_fasta, gff3_file,
                                   vcf_file, tmp_path):
        """Issue #2 edge case: gene with all-N ref completes with 0 codons."""
        result = runner.invoke(main, [
            "run", "--vcf", vcf_file, "--gff", gff3_file,
            "--fasta", ref_all_n_fasta, "--outdir", str(tmp_path),
            "--min-freq", "0", "--min-depth", "0", "--min-qual", "0",
        ])
        assert result.exit_code == 0, result.output
        df = pd.read_csv(tmp_path / "gene_results.tsv", sep="\t")
        g1 = df[df["gene_id"].str.contains("gene1", case=False)].iloc[0]
        assert g1["n_codons"] == 0

    def test_cdsonly_gff_exit_code_1(self, runner, ref_fasta, cdsonly_gff,
                                      vcf_file, tmp_path):
        """Issue #4: CDS-only GFF (no gene features) -> exit code 1."""
        result = runner.invoke(main, [
            "run", "--vcf", vcf_file, "--gff", cdsonly_gff,
            "--fasta", ref_fasta, "--outdir", str(tmp_path),
        ])
        assert result.exit_code == 1
        assert "No genes with CDS features found" in result.output

    def test_contig_mismatch_exit_code_1(self, runner, ref_fasta, gff3_file,
                                          mismatch_vcf_file, tmp_path):
        """Issue #7: chr1 (GFF) vs 1 (VCF) -> exit code 1 with helpful message."""
        result = runner.invoke(main, [
            "run", "--vcf", mismatch_vcf_file, "--gff", gff3_file,
            "--fasta", ref_fasta, "--outdir", str(tmp_path),
        ])
        assert result.exit_code == 1
        assert "No contig names shared" in result.output

    def test_summary_uses_cds_snp_variants(self, runner, ref_fasta, gff3_file,
                                            vcf_file, tmp_path):
        """Issue #6: summary column renamed from total_variants to cds_snp_variants."""
        result = runner.invoke(main, [
            "run", "--vcf", vcf_file, "--gff", gff3_file,
            "--fasta", ref_fasta, "--outdir", str(tmp_path),
            "--min-freq", "0", "--min-depth", "0", "--min-qual", "0",
        ])
        assert result.exit_code == 0, result.output
        df = pd.read_csv(tmp_path / "summary.tsv", sep="\t")
        assert "cds_snp_variants" in df.columns
        assert "total_variants" not in df.columns
