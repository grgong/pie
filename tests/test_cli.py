import pytest
from click.testing import CliRunner
from pie.cli import main


@pytest.fixture
def runner():
    return CliRunner()


class TestRunCommand:
    def test_help(self, runner):
        result = runner.invoke(main, ["run", "--help"])
        assert result.exit_code == 0
        assert "--vcf" in result.output
        assert "--gff" in result.output
        assert "--fasta" in result.output

    def test_missing_required(self, runner):
        result = runner.invoke(main, ["run"])
        assert result.exit_code != 0

    def test_full_run(self, runner, ref_fasta, gff3_file, vcf_file, tmp_path):
        result = runner.invoke(main, [
            "run",
            "--vcf", vcf_file,
            "--gff", gff3_file,
            "--fasta", ref_fasta,
            "--outdir", str(tmp_path),
            "--min-freq", "0",
            "--min-depth", "0",
            "--min-qual", "0",
        ])
        assert result.exit_code == 0, result.output
        assert (tmp_path / "gene_results.tsv").exists()
        assert (tmp_path / "window_results.tsv").exists()
        assert (tmp_path / "summary.tsv").exists()

    def test_help_shows_keep_multiallelic(self, runner):
        result = runner.invoke(main, ["run", "--help"])
        assert result.exit_code == 0
        assert "--keep-multiallelic" in result.output

    def test_run_with_threading(self, runner, ref_fasta, gff3_file, vcf_file, tmp_path):
        result = runner.invoke(main, [
            "run",
            "--vcf", vcf_file,
            "--gff", gff3_file,
            "--fasta", ref_fasta,
            "--outdir", str(tmp_path),
            "--min-freq", "0",
            "--min-depth", "0",
            "--min-qual", "0",
            "--threads", "2",
        ])
        assert result.exit_code == 0, result.output


class TestPlotCommand:
    def test_help(self, runner):
        result = runner.invoke(main, ["plot", "--help"])
        assert result.exit_code == 0
        assert "--gene-results" in result.output

    def test_plot_from_results(self, runner, ref_fasta, gff3_file, vcf_file, tmp_path):
        # First run the analysis
        runner.invoke(main, [
            "run", "--vcf", vcf_file, "--gff", gff3_file, "--fasta", ref_fasta,
            "--outdir", str(tmp_path), "--min-freq", "0", "--min-depth", "0", "--min-qual", "0",
        ])
        # Then plot
        result = runner.invoke(main, [
            "plot",
            "--gene-results", str(tmp_path / "gene_results.tsv"),
            "--output", str(tmp_path / "manhattan.png"),
        ])
        assert result.exit_code == 0, result.output
        assert (tmp_path / "manhattan.png").exists()


class TestSummaryCommand:
    def test_help(self, runner):
        result = runner.invoke(main, ["summary", "--help"])
        assert result.exit_code == 0

    def test_print_summary(self, runner, ref_fasta, gff3_file, vcf_file, tmp_path):
        runner.invoke(main, [
            "run", "--vcf", vcf_file, "--gff", gff3_file, "--fasta", ref_fasta,
            "--outdir", str(tmp_path), "--min-freq", "0", "--min-depth", "0", "--min-qual", "0",
        ])
        result = runner.invoke(main, ["summary", str(tmp_path / "summary.tsv")])
        assert result.exit_code == 0
        assert "genome_piN" in result.output
