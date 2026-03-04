import pandas as pd
import pytest
from pie.plot import manhattan_plot, scatter_plot, histogram_plot, boxplot_plot, sliding_window_plot


def _make_gene_tsv(tmp_path, df):
    """Helper: write a DataFrame to gene_results.tsv and return the path string."""
    tsv = tmp_path / "gene_results.tsv"
    df.to_csv(tsv, sep="\t", index=False)
    return str(tsv)


@pytest.fixture
def basic_gene_df():
    return pd.DataFrame({
        "chrom": ["chr1", "chr1", "chr2", "chr2"],
        "gene_id": ["g1", "g2", "g3", "g4"],
        "start": [100, 500, 100, 800],
        "end": [200, 600, 300, 900],
        "n_codons": [30, 40, 25, 50],
        "piN": [0.01, 0.02, 0.015, 0.03],
        "piS": [0.02, 0.017, 0.019, 0.015],
        "piN_piS": [0.5, 1.2, 0.8, 2.0],
    })


class TestManhattanPlot:
    def test_creates_png(self, tmp_path, basic_gene_df):
        tsv = _make_gene_tsv(tmp_path, basic_gene_df)
        out = str(tmp_path / "manhattan.png")
        manhattan_plot(tsv, out)
        assert (tmp_path / "manhattan.png").exists()
        assert (tmp_path / "manhattan.png").stat().st_size > 0

    def test_custom_dimensions(self, tmp_path, basic_gene_df):
        tsv = _make_gene_tsv(tmp_path, basic_gene_df)
        out = str(tmp_path / "manhattan.png")
        manhattan_plot(tsv, out, width=20, height=8)
        assert (tmp_path / "manhattan.png").exists()

    def test_log_scale(self, tmp_path, basic_gene_df):
        tsv = _make_gene_tsv(tmp_path, basic_gene_df)
        out = str(tmp_path / "manhattan.png")
        manhattan_plot(tsv, out, log_scale=True)
        assert (tmp_path / "manhattan.png").exists()

    def test_label_top(self, tmp_path, basic_gene_df):
        tsv = _make_gene_tsv(tmp_path, basic_gene_df)
        out = str(tmp_path / "manhattan.png")
        manhattan_plot(tsv, out, label_top=2)
        assert (tmp_path / "manhattan.png").exists()

    def test_highlight_genes(self, tmp_path, basic_gene_df):
        tsv = _make_gene_tsv(tmp_path, basic_gene_df)
        out = str(tmp_path / "manhattan.png")
        manhattan_plot(tsv, out, highlight_genes=["g1", "g3"])
        assert (tmp_path / "manhattan.png").exists()

    def test_numeric_chrom_ids(self, tmp_path):
        df = pd.DataFrame({
            "chrom": [1, 1, 2, 2],
            "gene_id": ["g1", "g2", "g3", "g4"],
            "start": [100, 500, 100, 800],
            "end": [200, 600, 300, 900],
            "n_codons": [30, 40, 25, 50],
            "piN": [0.01, 0.02, 0.015, 0.03],
            "piS": [0.02, 0.017, 0.019, 0.015],
            "piN_piS": [0.5, 1.2, 0.8, 2.0],
        })
        tsv = _make_gene_tsv(tmp_path, df)
        out = str(tmp_path / "manhattan.png")
        manhattan_plot(tsv, out)
        assert (tmp_path / "manhattan.png").exists()

    def test_handles_all_na(self, tmp_path):
        df = pd.DataFrame({
            "chrom": ["chr1"],
            "gene_id": ["g1"],
            "start": [100],
            "end": [200],
            "n_codons": [30],
            "piN": [0.0],
            "piS": [0.0],
            "piN_piS": [float("nan")],
        })
        tsv = _make_gene_tsv(tmp_path, df)
        out = str(tmp_path / "manhattan.png")
        manhattan_plot(tsv, out)
        assert (tmp_path / "manhattan.png").exists()

    def test_pdf_output(self, tmp_path, basic_gene_df):
        tsv = _make_gene_tsv(tmp_path, basic_gene_df)
        out = str(tmp_path / "manhattan.pdf")
        manhattan_plot(tsv, out)
        assert (tmp_path / "manhattan.pdf").exists()


class TestScatterPlot:
    def test_creates_png(self, tmp_path, basic_gene_df):
        tsv = _make_gene_tsv(tmp_path, basic_gene_df)
        out = str(tmp_path / "scatter.png")
        scatter_plot(tsv, out)
        assert (tmp_path / "scatter.png").exists()
        assert (tmp_path / "scatter.png").stat().st_size > 0

    def test_color_by_chrom(self, tmp_path, basic_gene_df):
        tsv = _make_gene_tsv(tmp_path, basic_gene_df)
        out = str(tmp_path / "scatter.png")
        scatter_plot(tsv, out, color_by_chrom=True)
        assert (tmp_path / "scatter.png").exists()

    def test_handles_all_na(self, tmp_path):
        df = pd.DataFrame({
            "chrom": ["chr1"],
            "gene_id": ["g1"],
            "start": [100],
            "end": [200],
            "n_codons": [30],
            "piN": [float("nan")],
            "piS": [float("nan")],
            "piN_piS": [float("nan")],
        })
        tsv = _make_gene_tsv(tmp_path, df)
        out = str(tmp_path / "scatter.png")
        scatter_plot(tsv, out)
        assert (tmp_path / "scatter.png").exists()


class TestHistogramPlot:
    def test_creates_png(self, tmp_path, basic_gene_df):
        tsv = _make_gene_tsv(tmp_path, basic_gene_df)
        out = str(tmp_path / "histogram.png")
        histogram_plot(tsv, out)
        assert (tmp_path / "histogram.png").exists()
        assert (tmp_path / "histogram.png").stat().st_size > 0

    def test_handles_all_na(self, tmp_path):
        df = pd.DataFrame({
            "chrom": ["chr1"],
            "gene_id": ["g1"],
            "start": [100],
            "end": [200],
            "n_codons": [30],
            "piN": [0.0],
            "piS": [0.0],
            "piN_piS": [float("nan")],
        })
        tsv = _make_gene_tsv(tmp_path, df)
        out = str(tmp_path / "histogram.png")
        histogram_plot(tsv, out)
        assert (tmp_path / "histogram.png").exists()


class TestBoxplotPlot:
    def test_creates_png(self, tmp_path, basic_gene_df):
        tsv = _make_gene_tsv(tmp_path, basic_gene_df)
        out = str(tmp_path / "boxplot.png")
        boxplot_plot(tsv, out)
        assert (tmp_path / "boxplot.png").exists()
        assert (tmp_path / "boxplot.png").stat().st_size > 0

    def test_single_chromosome(self, tmp_path):
        df = pd.DataFrame({
            "chrom": ["chr1", "chr1", "chr1"],
            "gene_id": ["g1", "g2", "g3"],
            "start": [100, 300, 500],
            "end": [200, 400, 600],
            "n_codons": [30, 40, 25],
            "piN": [0.01, 0.02, 0.015],
            "piS": [0.02, 0.017, 0.019],
            "piN_piS": [0.5, 1.2, 0.8],
        })
        tsv = _make_gene_tsv(tmp_path, df)
        out = str(tmp_path / "boxplot.png")
        boxplot_plot(tsv, out)
        assert (tmp_path / "boxplot.png").exists()

    def test_handles_all_na(self, tmp_path):
        df = pd.DataFrame({
            "chrom": ["chr1"],
            "gene_id": ["g1"],
            "start": [100],
            "end": [200],
            "n_codons": [30],
            "piN": [0.0],
            "piS": [0.0],
            "piN_piS": [float("nan")],
        })
        tsv = _make_gene_tsv(tmp_path, df)
        out = str(tmp_path / "boxplot.png")
        boxplot_plot(tsv, out)
        assert (tmp_path / "boxplot.png").exists()


class TestSlidingWindowPlot:
    def test_creates_png(self, tmp_path):
        df = pd.DataFrame({
            "chrom": ["chr1"] * 5 + ["chr2"] * 5,
            "win_start": list(range(0, 500, 100)) * 2,
            "win_end": list(range(100, 600, 100)) * 2,
            "gene_id": ["g1"] * 5 + ["g2"] * 5,
            "piN_piS": [0.5, 0.8, 1.0, 1.2, 0.9, 1.5, 1.1, 0.7, 0.6, 1.3],
        })
        tsv = tmp_path / "window_results.tsv"
        df.to_csv(tsv, sep="\t", index=False)
        out = str(tmp_path / "sliding_window.png")
        sliding_window_plot(str(tsv), out)
        assert (tmp_path / "sliding_window.png").exists()
        assert (tmp_path / "sliding_window.png").stat().st_size > 0

    def test_handles_empty(self, tmp_path):
        df = pd.DataFrame({
            "chrom": pd.Series([], dtype=str),
            "win_start": pd.Series([], dtype=int),
            "win_end": pd.Series([], dtype=int),
            "gene_id": pd.Series([], dtype=str),
            "piN_piS": pd.Series([], dtype=float),
        })
        tsv = tmp_path / "window_results.tsv"
        df.to_csv(tsv, sep="\t", index=False)
        out = str(tmp_path / "sliding_window.png")
        sliding_window_plot(str(tsv), out)
        assert (tmp_path / "sliding_window.png").exists()
