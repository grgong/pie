import pytest
import pandas as pd
from pie.io import write_gene_results, write_window_results, write_summary
from pie.diversity import GeneResult, CodonResult


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
