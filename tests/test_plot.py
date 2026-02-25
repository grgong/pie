import pytest
import pandas as pd
from pie.plot import manhattan_plot


class TestManhattanPlot:
    def test_creates_png(self, tmp_path):
        df = pd.DataFrame({
            "chrom": ["chr1", "chr1", "chr2", "chr2"],
            "gene_id": ["g1", "g2", "g3", "g4"],
            "start": [100, 500, 100, 800],
            "end": [200, 600, 300, 900],
            "piN_piS": [0.5, 1.2, 0.8, 2.0],
        })
        tsv = tmp_path / "gene_results.tsv"
        df.to_csv(tsv, sep="\t", index=False)
        out = tmp_path / "manhattan.png"
        manhattan_plot(str(tsv), str(out))
        assert out.exists()
        assert out.stat().st_size > 0

    def test_custom_dimensions(self, tmp_path):
        df = pd.DataFrame({
            "chrom": ["chr1"],
            "gene_id": ["g1"],
            "start": [100],
            "end": [200],
            "piN_piS": [1.0],
        })
        tsv = tmp_path / "gene_results.tsv"
        df.to_csv(tsv, sep="\t", index=False)
        out = tmp_path / "manhattan.png"
        manhattan_plot(str(tsv), str(out), width=20, height=8)
        assert out.exists()

    def test_handles_all_na(self, tmp_path):
        df = pd.DataFrame({
            "chrom": ["chr1"],
            "gene_id": ["g1"],
            "start": [100],
            "end": [200],
            "piN_piS": [float("nan")],
        })
        tsv = tmp_path / "gene_results.tsv"
        df.to_csv(tsv, sep="\t", index=False)
        out = tmp_path / "manhattan.png"
        manhattan_plot(str(tsv), str(out))
        assert out.exists()
