from pie.cli import main


class TestRunCommand:
    def test_run_group_help(self, runner):
        result = runner.invoke(main, ["run", "--help"])
        assert result.exit_code == 0
        assert "pool" in result.output
        assert "ind" in result.output

    def test_run_no_subcommand_shows_help(self, runner):
        result = runner.invoke(main, ["run"])
        assert "pool" in result.output
        assert "ind" in result.output

    def test_pool_help(self, runner):
        result = runner.invoke(main, ["run", "pool", "--help"])
        assert result.exit_code == 0
        assert "--vcf" in result.output
        assert "--gff" in result.output
        assert "--fasta" in result.output
        assert "--min-depth" in result.output
        assert "--sample" in result.output

    def test_ind_help(self, runner):
        result = runner.invoke(main, ["run", "ind", "--help"])
        assert result.exit_code == 0
        assert "--vcf" in result.output
        assert "--samples" in result.output
        assert "--samples-file" in result.output
        assert "--min-call-rate" in result.output
        assert "--min-an" in result.output

    def test_pool_missing_required(self, runner):
        result = runner.invoke(main, ["run", "pool"])
        assert result.exit_code != 0

    def test_full_run_pool(self, runner, ref_fasta, gff3_file, vcf_file, tmp_path):
        result = runner.invoke(main, [
            "run", "pool",
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

    def test_pool_help_shows_keep_multiallelic(self, runner):
        result = runner.invoke(main, ["run", "pool", "--help"])
        assert result.exit_code == 0
        assert "--keep-multiallelic" in result.output

    def test_pool_help_shows_quiet(self, runner):
        result = runner.invoke(main, ["run", "pool", "--help"])
        assert result.exit_code == 0
        assert "--quiet" in result.output

    def test_quiet_flag_runs(self, runner, ref_fasta, gff3_file, vcf_file, tmp_path):
        """--quiet flag is accepted and suppresses INFO messages."""
        result = runner.invoke(main, [
            "run", "pool",
            "--vcf", vcf_file,
            "--gff", gff3_file,
            "--fasta", ref_fasta,
            "--outdir", str(tmp_path),
            "--min-freq", "0",
            "--min-depth", "0",
            "--min-qual", "0",
            "--quiet",
        ])
        assert result.exit_code == 0, result.output
        assert (tmp_path / "gene_results.tsv").exists()

    def test_pool_help_shows_include_stop_codons(self, runner):
        result = runner.invoke(main, ["run", "pool", "--help"])
        assert result.exit_code == 0
        assert "--include-stop-codons" in result.output

    def test_include_stop_codons_flag(self, runner, ref_fasta, gff3_file, vcf_file, tmp_path):
        """--include-stop-codons flag is accepted and runs successfully."""
        result = runner.invoke(main, [
            "run", "pool",
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
            "run", "pool",
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
    def test_plot_group_help(self, runner):
        result = runner.invoke(main, ["plot", "--help"])
        assert result.exit_code == 0
        assert "manhattan" in result.output
        assert "scatter" in result.output
        assert "histogram" in result.output
        assert "boxplot" in result.output
        assert "sliding-window" in result.output

    def test_manhattan_help(self, runner):
        result = runner.invoke(main, ["plot", "manhattan", "--help"])
        assert result.exit_code == 0
        assert "--input" in result.output
        assert "--output" in result.output
        assert "--log-scale" in result.output
        assert "--label-top" in result.output
        assert "--highlight-genes" in result.output

    def test_scatter_help(self, runner):
        result = runner.invoke(main, ["plot", "scatter", "--help"])
        assert result.exit_code == 0
        assert "--color-by-chrom" in result.output

    def _run_pool(self, runner, vcf_file, gff3_file, ref_fasta, tmp_path):
        runner.invoke(main, [
            "run", "pool", "--vcf", vcf_file, "--gff", gff3_file, "--fasta", ref_fasta,
            "--outdir", str(tmp_path), "--min-freq", "0", "--min-depth", "0", "--min-qual", "0",
        ])

    def test_manhattan_from_results(self, runner, ref_fasta, gff3_file, vcf_file, tmp_path):
        self._run_pool(runner, vcf_file, gff3_file, ref_fasta, tmp_path)
        result = runner.invoke(main, [
            "plot", "manhattan",
            "-i", str(tmp_path / "gene_results.tsv"),
            "-o", str(tmp_path / "manhattan.png"),
        ])
        assert result.exit_code == 0, result.output
        assert (tmp_path / "manhattan.png").exists()

    def test_scatter_from_results(self, runner, ref_fasta, gff3_file, vcf_file, tmp_path):
        self._run_pool(runner, vcf_file, gff3_file, ref_fasta, tmp_path)
        result = runner.invoke(main, [
            "plot", "scatter",
            "-i", str(tmp_path / "gene_results.tsv"),
            "-o", str(tmp_path / "scatter.png"),
        ])
        assert result.exit_code == 0, result.output
        assert (tmp_path / "scatter.png").exists()

    def test_boxplot_from_results(self, runner, ref_fasta, gff3_file, vcf_file, tmp_path):
        self._run_pool(runner, vcf_file, gff3_file, ref_fasta, tmp_path)
        result = runner.invoke(main, [
            "plot", "boxplot",
            "-i", str(tmp_path / "gene_results.tsv"),
            "-o", str(tmp_path / "boxplot.png"),
        ])
        assert result.exit_code == 0, result.output
        assert (tmp_path / "boxplot.png").exists()

    def test_histogram_from_results(self, runner, ref_fasta, gff3_file, vcf_file, tmp_path):
        self._run_pool(runner, vcf_file, gff3_file, ref_fasta, tmp_path)
        result = runner.invoke(main, [
            "plot", "histogram",
            "-i", str(tmp_path / "gene_results.tsv"),
            "-o", str(tmp_path / "histogram.png"),
        ])
        assert result.exit_code == 0, result.output
        assert (tmp_path / "histogram.png").exists()

    def test_sliding_window_from_results(self, runner, ref_fasta, gff3_file, vcf_file, tmp_path):
        self._run_pool(runner, vcf_file, gff3_file, ref_fasta, tmp_path)
        result = runner.invoke(main, [
            "plot", "sliding-window",
            "-i", str(tmp_path / "window_results.tsv"),
            "-o", str(tmp_path / "sw.png"),
        ])
        assert result.exit_code == 0, result.output
        assert (tmp_path / "sw.png").exists()

    def test_filter_options_in_help(self, runner):
        """All plot subcommands expose shared and metric-specific filter options."""
        for cmd, expected in [
            ("manhattan", ["--max-ratio", "--exclude-zero-ratio", "--min-codons", "--min-variants"]),
            ("scatter", ["--max-piN", "--max-piS", "--min-codons", "--min-variants"]),
            ("histogram", ["--max-ratio", "--exclude-zero-ratio", "--min-codons", "--min-variants"]),
            ("boxplot", ["--max-piN", "--max-piS", "--max-ratio", "--exclude-zero-ratio", "--min-codons", "--min-variants"]),
            ("sliding-window", ["--max-ratio"]),
        ]:
            result = runner.invoke(main, ["plot", cmd, "--help"])
            assert result.exit_code == 0
            for opt in expected:
                assert opt in result.output, f"{opt} missing from 'plot {cmd} --help'"

    def test_manhattan_with_filters(self, runner, ref_fasta, gff3_file, vcf_file, tmp_path):
        self._run_pool(runner, vcf_file, gff3_file, ref_fasta, tmp_path)
        result = runner.invoke(main, [
            "plot", "manhattan",
            "-i", str(tmp_path / "gene_results.tsv"),
            "-o", str(tmp_path / "manhattan.png"),
            "--max-ratio", "1.5", "--exclude-zero-ratio",
        ])
        assert result.exit_code == 0, result.output

    def test_boxplot_with_filters(self, runner, ref_fasta, gff3_file, vcf_file, tmp_path):
        self._run_pool(runner, vcf_file, gff3_file, ref_fasta, tmp_path)
        result = runner.invoke(main, [
            "plot", "boxplot",
            "-i", str(tmp_path / "gene_results.tsv"),
            "-o", str(tmp_path / "boxplot.png"),
            "--max-piN", "1.0", "--max-piS", "1.0", "--max-ratio", "1.5",
        ])
        assert result.exit_code == 0, result.output


class TestSummaryCommand:
    def test_help(self, runner):
        result = runner.invoke(main, ["summary", "--help"])
        assert result.exit_code == 0

    def test_print_summary(self, runner, ref_fasta, gff3_file, vcf_file, tmp_path):
        runner.invoke(main, [
            "run", "pool", "--vcf", vcf_file, "--gff", gff3_file, "--fasta", ref_fasta,
            "--outdir", str(tmp_path), "--min-freq", "0", "--min-depth", "0", "--min-qual", "0",
        ])
        result = runner.invoke(main, ["summary", str(tmp_path / "summary.tsv")])
        assert result.exit_code == 0
        assert "genome_piN" in result.output


class TestIndividualMode:
    def test_ind_runs(self, runner, ref_fasta, gff3_file, individual_vcf_file, tmp_path):
        """pie run ind works end-to-end."""
        result = runner.invoke(main, [
            "run", "ind",
            "--vcf", individual_vcf_file,
            "--gff", gff3_file, "--fasta", ref_fasta,
            "--outdir", str(tmp_path),
            "--min-freq", "0", "--min-qual", "0", "--min-call-rate", "0",
        ])
        assert result.exit_code == 0, result.output

    def test_samples_and_samples_file_mutually_exclusive(self, runner, ref_fasta, gff3_file,
                                                          individual_vcf_file, tmp_path):
        sf = tmp_path / "samples.txt"
        sf.write_text("S1\nS2\n")
        result = runner.invoke(main, [
            "run", "ind",
            "--samples", "S1,S2", "--samples-file", str(sf),
            "--vcf", individual_vcf_file, "--gff", gff3_file,
            "--fasta", ref_fasta, "--outdir", str(tmp_path),
        ])
        assert result.exit_code != 0
        assert "mutually exclusive" in result.output.lower()

    def test_min_call_rate_out_of_range(self, runner, ref_fasta, gff3_file, vcf_file, tmp_path):
        result = runner.invoke(main, [
            "run", "ind", "--min-call-rate", "1.5",
            "--vcf", vcf_file, "--gff", gff3_file,
            "--fasta", ref_fasta, "--outdir", str(tmp_path),
        ])
        assert result.exit_code != 0

    def test_individual_no_samples_uses_all(self, runner, ref_fasta, gff3_file,
                                             individual_vcf_file, tmp_path):
        """pie run ind without --samples -> uses all samples in VCF."""
        result = runner.invoke(main, [
            "run", "ind",
            "--vcf", individual_vcf_file, "--gff", gff3_file,
            "--fasta", ref_fasta, "--outdir", str(tmp_path),
            "--min-freq", "0", "--min-qual", "0", "--min-call-rate", "0",
        ])
        assert result.exit_code == 0, result.output

    def test_invalid_sample_name_error(self, runner, ref_fasta, gff3_file,
                                       individual_vcf_file, tmp_path):
        """--samples with name not in VCF -> error with available names."""
        result = runner.invoke(main, [
            "run", "ind",
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
            "run", "ind",
            "--samples-file", str(sf),
            "--vcf", individual_vcf_file, "--gff", gff3_file,
            "--fasta", ref_fasta, "--outdir", str(tmp_path),
            "--min-freq", "0", "--min-qual", "0", "--min-call-rate", "0",
        ])
        assert result.exit_code == 0, result.output
