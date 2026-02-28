from pie.cli import main


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

    def test_help_shows_include_stop_codons(self, runner):
        result = runner.invoke(main, ["run", "--help"])
        assert result.exit_code == 0
        assert "--include-stop-codons" in result.output

    def test_include_stop_codons_flag(self, runner, ref_fasta, gff3_file, vcf_file, tmp_path):
        """--include-stop-codons flag is accepted and runs successfully."""
        result = runner.invoke(main, [
            "run",
            "--vcf", vcf_file,
            "--gff", gff3_file,
            "--fasta", ref_fasta,
            "--outdir", str(tmp_path),
            "--min-freq", "0",
            "--min-depth", "0",
            "--min-qual", "0",
            "--include-stop-codons",
        ])
        assert result.exit_code == 0, result.output
        assert (tmp_path / "gene_results.tsv").exists()

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


class TestModeValidation:
    def test_help_shows_mode(self, runner):
        result = runner.invoke(main, ["run", "--help"])
        assert "--mode" in result.output

    def test_default_mode_is_pool(self, runner, ref_fasta, gff3_file, vcf_file, tmp_path):
        """No --mode -> defaults to pool, works as before."""
        result = runner.invoke(main, [
            "run", "--vcf", vcf_file, "--gff", gff3_file,
            "--fasta", ref_fasta, "--outdir", str(tmp_path),
            "--min-freq", "0", "--min-depth", "0", "--min-qual", "0",
        ])
        assert result.exit_code == 0, result.output

    def test_mode_pool_explicit(self, runner, ref_fasta, gff3_file, vcf_file, tmp_path):
        result = runner.invoke(main, [
            "run", "--mode", "pool", "--vcf", vcf_file, "--gff", gff3_file,
            "--fasta", ref_fasta, "--outdir", str(tmp_path),
            "--min-freq", "0", "--min-depth", "0", "--min-qual", "0",
        ])
        assert result.exit_code == 0, result.output

    def test_mode_ind_alias(self, runner, ref_fasta, gff3_file, individual_vcf_file, tmp_path):
        """'ind' is accepted as alias for 'individual'."""
        result = runner.invoke(main, [
            "run", "--mode", "ind", "--vcf", individual_vcf_file,
            "--gff", gff3_file, "--fasta", ref_fasta,
            "--outdir", str(tmp_path),
            "--min-freq", "0", "--min-qual", "0", "--min-call-rate", "0",
        ])
        assert result.exit_code == 0, result.output

    def test_pool_with_samples_error(self, runner, ref_fasta, gff3_file, vcf_file, tmp_path):
        """--mode pool + --samples -> error."""
        result = runner.invoke(main, [
            "run", "--mode", "pool", "--samples", "S1,S2",
            "--vcf", vcf_file, "--gff", gff3_file,
            "--fasta", ref_fasta, "--outdir", str(tmp_path),
        ])
        assert result.exit_code != 0
        assert "individual" in result.output.lower() or "pool" in result.output.lower()

    def test_pool_with_samples_file_error(self, runner, ref_fasta, gff3_file, vcf_file, tmp_path):
        """--mode pool + --samples-file -> error."""
        sf = tmp_path / "samples.txt"
        sf.write_text("S1\nS2\n")
        result = runner.invoke(main, [
            "run", "--mode", "pool", "--samples-file", str(sf),
            "--vcf", vcf_file, "--gff", gff3_file,
            "--fasta", ref_fasta, "--outdir", str(tmp_path),
        ])
        assert result.exit_code != 0

    def test_pool_with_min_call_rate_error(self, runner, ref_fasta, gff3_file, vcf_file, tmp_path):
        """--mode pool + --min-call-rate -> error."""
        result = runner.invoke(main, [
            "run", "--mode", "pool", "--min-call-rate", "0.5",
            "--vcf", vcf_file, "--gff", gff3_file,
            "--fasta", ref_fasta, "--outdir", str(tmp_path),
        ])
        assert result.exit_code != 0

    def test_pool_with_min_an_error(self, runner, ref_fasta, gff3_file, vcf_file, tmp_path):
        """--mode pool + --min-an -> error."""
        result = runner.invoke(main, [
            "run", "--mode", "pool", "--min-an", "4",
            "--vcf", vcf_file, "--gff", gff3_file,
            "--fasta", ref_fasta, "--outdir", str(tmp_path),
        ])
        assert result.exit_code != 0

    def test_individual_with_sample_error(self, runner, ref_fasta, gff3_file, vcf_file, tmp_path):
        """--mode individual + --sample -> error."""
        result = runner.invoke(main, [
            "run", "--mode", "individual", "--sample", "S1",
            "--vcf", vcf_file, "--gff", gff3_file,
            "--fasta", ref_fasta, "--outdir", str(tmp_path),
        ])
        assert result.exit_code != 0

    def test_individual_with_min_depth_error(self, runner, ref_fasta, gff3_file, vcf_file, tmp_path):
        """--mode individual + --min-depth -> error."""
        result = runner.invoke(main, [
            "run", "--mode", "individual", "--min-depth", "10",
            "--vcf", vcf_file, "--gff", gff3_file,
            "--fasta", ref_fasta, "--outdir", str(tmp_path),
        ])
        assert result.exit_code != 0

    def test_samples_and_samples_file_mutually_exclusive(self, runner, ref_fasta, gff3_file,
                                                          individual_vcf_file, tmp_path):
        sf = tmp_path / "samples.txt"
        sf.write_text("S1\nS2\n")
        result = runner.invoke(main, [
            "run", "--mode", "individual",
            "--samples", "S1,S2", "--samples-file", str(sf),
            "--vcf", individual_vcf_file, "--gff", gff3_file,
            "--fasta", ref_fasta, "--outdir", str(tmp_path),
        ])
        assert result.exit_code != 0
        assert "mutually exclusive" in result.output.lower()

    def test_min_call_rate_out_of_range(self, runner, ref_fasta, gff3_file, vcf_file, tmp_path):
        result = runner.invoke(main, [
            "run", "--mode", "individual", "--min-call-rate", "1.5",
            "--vcf", vcf_file, "--gff", gff3_file,
            "--fasta", ref_fasta, "--outdir", str(tmp_path),
        ])
        assert result.exit_code != 0

    def test_individual_no_samples_uses_all(self, runner, ref_fasta, gff3_file,
                                             individual_vcf_file, tmp_path):
        """--mode individual without --samples -> uses all samples in VCF."""
        result = runner.invoke(main, [
            "run", "--mode", "individual",
            "--vcf", individual_vcf_file, "--gff", gff3_file,
            "--fasta", ref_fasta, "--outdir", str(tmp_path),
            "--min-freq", "0", "--min-qual", "0", "--min-call-rate", "0",
        ])
        assert result.exit_code == 0, result.output

    def test_invalid_sample_name_error(self, runner, ref_fasta, gff3_file,
                                       individual_vcf_file, tmp_path):
        """--samples with name not in VCF -> error with available names."""
        result = runner.invoke(main, [
            "run", "--mode", "individual",
            "--samples", "S1,NONEXISTENT",
            "--vcf", individual_vcf_file, "--gff", gff3_file,
            "--fasta", ref_fasta, "--outdir", str(tmp_path),
        ])
        assert result.exit_code != 0
        assert "NONEXISTENT" in result.output

    def test_samples_file(self, runner, ref_fasta, gff3_file,
                          individual_vcf_file, tmp_path):
        """--samples-file reads sample names from file."""
        sf = tmp_path / "samples.txt"
        sf.write_text("S1\nS2\n")
        result = runner.invoke(main, [
            "run", "--mode", "individual",
            "--samples-file", str(sf),
            "--vcf", individual_vcf_file, "--gff", gff3_file,
            "--fasta", ref_fasta, "--outdir", str(tmp_path),
            "--min-freq", "0", "--min-qual", "0", "--min-call-rate", "0",
        ])
        assert result.exit_code == 0, result.output
