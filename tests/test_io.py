import pytest
import pandas as pd
from pie.io import write_gene_results, write_window_results, write_summary, write_variant_results
from pie.diversity import GeneResult, CodonResult, VariantRecord


class TestIndividualModeColumns:
    def _make_result(self, n_samples=None, call_rates=None):
        """Helper to create a GeneResult with optional ind-mode metadata."""
        return GeneResult(
            gene_id="gene1", transcript_id="tx1", chrom="chr1",
            start=0, end=90, strand="+", n_codons=29, n_poly_codons=2,
            N_sites=59.0, S_sites=27.0, N_diffs=0.42, S_diffs=0.32,
            mean_variant_depth=100.0, n_variants=2,
            n_samples=n_samples, call_rates=call_rates,
        )

    def test_pool_mode_no_extra_columns(self, tmp_path):
        """Pool mode: n_samples=None -> no n_samples/mean_call_rate columns."""
        results = [self._make_result()]
        path = str(tmp_path / "gene_results.tsv")
        write_gene_results(results, path)
        df = pd.read_csv(path, sep="\t")
        assert "n_samples" not in df.columns
        assert "mean_call_rate" not in df.columns

    def test_individual_mode_extra_columns(self, tmp_path):
        """Individual mode: n_samples set -> extra columns present."""
        results = [self._make_result(n_samples=4, call_rates=[0.75, 1.0])]
        path = str(tmp_path / "gene_results.tsv")
        write_gene_results(results, path)
        df = pd.read_csv(path, sep="\t")
        assert "n_samples" in df.columns
        assert "mean_call_rate" in df.columns
        assert df.iloc[0]["n_samples"] == 4
        assert abs(df.iloc[0]["mean_call_rate"] - 0.875) < 1e-6  # (0.75+1.0)/2

    def test_summary_pool_mode(self, tmp_path):
        """Pool mode summary: no n_samples_selected/mean_call_rate."""
        results = [self._make_result()]
        path = str(tmp_path / "summary.tsv")
        write_summary(results, path)
        df = pd.read_csv(path, sep="\t")
        assert "n_samples_selected" not in df.columns
        assert "mean_call_rate" not in df.columns

    def test_summary_individual_mode(self, tmp_path):
        """Individual mode summary: includes n_samples_selected and mean_call_rate."""
        r1 = self._make_result(n_samples=4, call_rates=[0.75, 1.0])
        r2 = self._make_result(n_samples=4, call_rates=[1.0])
        results = [r1, r2]
        path = str(tmp_path / "summary.tsv")
        write_summary(results, path)
        df = pd.read_csv(path, sep="\t")
        assert df.iloc[0]["n_samples_selected"] == 4
        # Variant-site-weighted: (0.75 + 1.0 + 1.0) / 3 = 0.9167
        assert abs(df.iloc[0]["mean_call_rate"] - (0.75 + 1.0 + 1.0) / 3) < 1e-4


@pytest.fixture
def sample_results():
    return [
        GeneResult(
            gene_id="gene1", transcript_id="tx1", chrom="chr1",
            start=0, end=90, strand="+", n_codons=29, n_poly_codons=2,
            N_sites=59.6667, S_sites=27.3333,
            N_diffs=0.42, S_diffs=0.32,
            mean_variant_depth=100.0, n_variants=2,
            codon_results=[
                CodonResult("chr1", 0, 2.5, 0.5, 0.0, 0.0),
                CodonResult("chr1", 3, 2.0, 1.0, 0.0, 0.32),
                CodonResult("chr1", 6, 8/3, 1/3, 0.42, 0.0),
                CodonResult("chr1", 9, 2.0, 1.0, 0.0, 0.0),
            ],
        ),
        GeneResult(
            gene_id="gene2", transcript_id="tx2", chrom="chr1",
            start=100, end=220, strand="+", n_codons=32, n_poly_codons=1,
            N_sites=65.6667, S_sites=30.3333,
            N_diffs=0.48, S_diffs=0.0,
            mean_variant_depth=100.0, n_variants=1,
            codon_results=[CodonResult("chr1", 100, 2.0, 1.0, 0.0, 0.0)],
        ),
    ]


class TestWriteGeneResults:
    def test_writes_tsv(self, sample_results, tmp_path):
        outpath = tmp_path / "gene_results.tsv"
        write_gene_results(sample_results, str(outpath))
        assert outpath.exists()
        df = pd.read_csv(outpath, sep="\t")
        assert len(df) == 2
        assert set(df.columns) >= {"piN", "piS", "piN_piS", "gene_id"}

    def test_correct_piN(self, sample_results, tmp_path):
        outpath = tmp_path / "gene_results.tsv"
        write_gene_results(sample_results, str(outpath))
        df = pd.read_csv(outpath, sep="\t")
        row = df[df.gene_id == "gene1"].iloc[0]
        assert abs(row["piN"] - 0.42/59.6667) < 1e-6

    def test_piN_piS_na(self, sample_results, tmp_path):
        outpath = tmp_path / "gene_results.tsv"
        write_gene_results(sample_results, str(outpath))
        df = pd.read_csv(outpath, sep="\t")
        row = df[df.gene_id == "gene2"].iloc[0]
        assert pd.isna(row["piN_piS"])  # piS = 0


class TestWriteWindowResults:
    def test_writes_tsv(self, sample_results, tmp_path):
        outpath = tmp_path / "window_results.tsv"
        write_window_results(sample_results, str(outpath), window_size=10, window_step=5)
        assert outpath.exists()
        df = pd.read_csv(outpath, sep="\t")
        assert len(df) > 0
        assert "piN" in df.columns


class TestWriteSummary:
    def test_writes_tsv(self, sample_results, tmp_path):
        outpath = tmp_path / "summary.tsv"
        write_summary(sample_results, str(outpath))
        assert outpath.exists()
        df = pd.read_csv(outpath, sep="\t")
        assert len(df) == 1
        assert "genome_piN" in df.columns
        # genome_piN = (0.42 + 0.48) / (59.6667 + 65.6667)
        expected_piN = 0.9 / 125.3334
        assert abs(df.iloc[0]["genome_piN"] - expected_piN) < 1e-6

    def test_stop_renorm_columns_no_stops(self, sample_results, tmp_path):
        """No stop-codon renormalization -> both columns are 0."""
        outpath = tmp_path / "summary.tsv"
        write_summary(sample_results, str(outpath))
        df = pd.read_csv(outpath, sep="\t")
        assert df.iloc[0]["stop_renorm_genes"] == 0
        assert df.iloc[0]["stop_renorm_codons"] == 0

    def test_stop_renorm_columns_with_stops(self, tmp_path):
        """Genes with n_stop_codons > 0 are counted correctly."""
        r1 = GeneResult(
            gene_id="g1", transcript_id="tx1", chrom="chr1",
            start=0, end=90, strand="+", n_codons=30, n_poly_codons=2,
            N_sites=60.0, S_sites=27.0, N_diffs=0.4, S_diffs=0.3,
            mean_variant_depth=100.0, n_variants=2, n_stop_codons=3,
        )
        r2 = GeneResult(
            gene_id="g2", transcript_id="tx2", chrom="chr1",
            start=100, end=190, strand="+", n_codons=30, n_poly_codons=1,
            N_sites=60.0, S_sites=27.0, N_diffs=0.5, S_diffs=0.0,
            mean_variant_depth=100.0, n_variants=1, n_stop_codons=0,
        )
        r3 = GeneResult(
            gene_id="g3", transcript_id="tx3", chrom="chr1",
            start=200, end=290, strand="+", n_codons=30, n_poly_codons=1,
            N_sites=60.0, S_sites=27.0, N_diffs=0.2, S_diffs=0.1,
            mean_variant_depth=100.0, n_variants=1, n_stop_codons=1,
        )
        outpath = tmp_path / "summary.tsv"
        write_summary([r1, r2, r3], str(outpath))
        df = pd.read_csv(outpath, sep="\t")
        assert df.iloc[0]["stop_renorm_genes"] == 2   # g1 and g3
        assert df.iloc[0]["stop_renorm_codons"] == 4   # 3 + 0 + 1


class TestWriteVariantResults:
    def _make_records(self):
        return [
            VariantRecord(
                chrom="chr1", pos=6, ref="T", alt="C", gene_id="gene1",
                codon_pos=3, ref_codon="GCT", alt_codon="GCC",
                ref_aa="A", alt_aa="A", variant_class="synonymous",
                ao=20, ro=80, dp=100, af=0.2,
                strand="+", cds_position=6, n_sites=0.0, s_sites=1.0,
            ),
            VariantRecord(
                chrom="chr1", pos=7, ref="G", alt="A", gene_id="gene1",
                codon_pos=1, ref_codon="GAT", alt_codon="AAT",
                ref_aa="D", alt_aa="N", variant_class="nonsynonymous",
                ao=30, ro=70, dp=100, af=0.3,
                strand="+", cds_position=7, n_sites=1.0, s_sites=0.0,
            ),
        ]

    def test_writes_tsv(self, tmp_path):
        path = str(tmp_path / "variant_results.tsv")
        write_variant_results(self._make_records(), path)
        df = pd.read_csv(path, sep="\t")
        assert len(df) == 2
        assert set(df.columns) >= {
            "chrom", "pos", "ref", "alt", "gene_id", "variant_class",
            "ao", "ro", "dp", "af",
        }

    def test_correct_values(self, tmp_path):
        path = str(tmp_path / "variant_results.tsv")
        write_variant_results(self._make_records(), path)
        df = pd.read_csv(path, sep="\t")
        row = df.iloc[0]
        assert row["chrom"] == "chr1"
        assert row["pos"] == 6
        assert row["ref"] == "T"
        assert row["alt"] == "C"
        assert row["variant_class"] == "synonymous"
        assert row["ao"] == 20
        assert row["ro"] == 80
        assert abs(row["af"] - 0.2) < 1e-6

    def test_empty_records(self, tmp_path):
        path = str(tmp_path / "variant_results.tsv")
        write_variant_results([], path)
        df = pd.read_csv(path, sep="\t")
        assert len(df) == 0
        assert "chrom" in df.columns


class TestEmptyWindowsSkipped:
    """Windows step on genomic position, codons are selected by position, so a
    window can land entirely inside an intron. Writing those as piN=piS=0.0 is
    indistinguishable from a real measurement of zero — and at the default
    100 bp step a gene with a 100 kb intron emits ~1000 of them."""

    @pytest.fixture
    def intron_gene(self):
        """10 codons at 1000-1029, a 20 kb intron, 10 codons at 21000-21029."""
        codons = [CodonResult("chr1", p, 2.0, 1.0, 0.0, 0.0)
                  for p in list(range(1000, 1030, 3)) + list(range(21000, 21030, 3))]
        return [GeneResult(
            gene_id="gene1", transcript_id="tx1", chrom="chr1",
            start=1000, end=21030, strand="+", n_codons=20, n_poly_codons=0,
            N_sites=40.0, S_sites=20.0, N_diffs=0.0, S_diffs=0.0,
            mean_variant_depth=0.0, n_variants=0, codon_results=codons,
        )]

    @pytest.fixture
    def window_df(self, intron_gene, tmp_path):
        outpath = tmp_path / "window_results.tsv"
        write_window_results(intron_gene, str(outpath),
                             window_size=1000, window_step=1000)
        return pd.read_csv(outpath, sep="\t")

    def test_intron_windows_are_not_written(self, window_df):
        assert len(window_df) == 2  # one per CDS block, not 21
        assert (window_df["n_codons"] > 0).all()

    def test_windows_keep_a_uniform_width(self, window_df):
        """plot.py derives win_mid from the interval; clipping would skew it."""
        assert (window_df["win_end"] - window_df["win_start"] == 1000).all()


class TestEmptyResultsKeepHeader:
    """pd.DataFrame([]) has no columns, so to_csv writes a bare newline."""

    def test_empty_window_results(self, tmp_path):
        outpath = tmp_path / "window_results.tsv"
        write_window_results([], str(outpath), window_size=1000, window_step=100)
        df = pd.read_csv(outpath, sep="\t")
        assert len(df) == 0
        assert {"chrom", "win_start", "win_end", "piN", "piS"} <= set(df.columns)

    def test_empty_gene_results(self, tmp_path):
        outpath = tmp_path / "gene_results.tsv"
        write_gene_results([], str(outpath))
        df = pd.read_csv(outpath, sep="\t")
        assert len(df) == 0
        assert {"gene_id", "piN", "piS", "piN_piS"} <= set(df.columns)

    def test_empty_header_matches_populated_header(self, sample_results,
                                                   tmp_path):
        """The declared empty-case columns must not drift from the row keys."""
        empty, populated = tmp_path / "e.tsv", tmp_path / "p.tsv"
        write_gene_results([], str(empty))
        write_gene_results(sample_results, str(populated))
        assert (list(pd.read_csv(empty, sep="\t").columns)
                == list(pd.read_csv(populated, sep="\t").columns))

    def test_empty_window_header_matches_populated(self, sample_results,
                                                   tmp_path):
        empty, populated = tmp_path / "e.tsv", tmp_path / "p.tsv"
        write_window_results([], str(empty), 1000, 100)
        write_window_results(sample_results, str(populated), 1000, 100)
        assert (list(pd.read_csv(empty, sep="\t").columns)
                == list(pd.read_csv(populated, sep="\t").columns))

    def test_gene_column_order_is_stable(self, sample_results, tmp_path):
        outpath = tmp_path / "gene_results.tsv"
        write_gene_results(sample_results, str(outpath))
        df = pd.read_csv(outpath, sep="\t")
        assert list(df.columns)[:6] == [
            "chrom", "gene_id", "transcript_id", "start", "end", "strand"]
