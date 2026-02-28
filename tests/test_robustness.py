"""Regression tests for robustness branches added in commit 0c878fd.

Covers:
  1. NoGenesFoundError when GFF has no genes with CDS
  2. Zero shared contigs → ValueError fail-fast
  3. Partial contig mismatch → genes on missing contigs dropped
  4. Ambiguous reference codon (N bases) filtering
  5. All-ambiguous codons → empty GeneResult
  6. summary.tsv uses cds_snp_variants (not total_variants)
"""

import subprocess
import textwrap

import pandas as pd
import pytest
from click.testing import CliRunner

from pie.annotation import NoGenesFoundError
from pie.cli import main
from pie.diversity import GeneResult, compute_gene_diversity
from pie.parallel import run_parallel


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _bgzip_and_index(vcf_path):
    """Bgzip and tabix-index a plain VCF file.  Returns .vcf.gz path."""
    gz_path = str(vcf_path) + ".gz"
    with open(gz_path, "wb") as out:
        subprocess.run(["bgzip", "-c", str(vcf_path)], stdout=out, check=True)
    subprocess.run(["tabix", "-p", "vcf", gz_path], check=True)
    return gz_path


def _write_fasta(path, sequences: dict[str, str]):
    """Write a multi-sequence FASTA and samtools-index it."""
    with open(path, "w") as fh:
        for name, seq in sequences.items():
            fh.write(f">{name}\n{seq}\n")
    subprocess.run(["samtools", "faidx", str(path)], check=True)
    return str(path)


@pytest.fixture
def runner():
    return CliRunner()


# ---------------------------------------------------------------------------
# 1. NoGenesFoundError — GFF with no CDS features
# ---------------------------------------------------------------------------

class TestNoGenesFoundError:
    def test_gff_without_cds_raises_error(self, runner, ref_fasta, vcf_file,
                                           tmp_path):
        """GFF containing only gene/mRNA but no CDS → exit 1 with message."""
        gff_no_cds = tmp_path / "no_cds.gff3"
        gff_no_cds.write_text(textwrap.dedent("""\
            ##gff-version 3
            chr1\ttest\tgene\t1\t90\t.\t+\t.\tID=gene1
            chr1\ttest\tmRNA\t1\t90\t.\t+\t.\tID=mRNA1;Parent=gene1
        """))
        result = runner.invoke(main, [
            "run", "--vcf", vcf_file, "--gff", str(gff_no_cds),
            "--fasta", ref_fasta, "--outdir", str(tmp_path / "out"),
            "--min-freq", "0", "--min-depth", "0", "--min-qual", "0",
        ])
        assert result.exit_code == 1
        assert "No genes" in result.output or "CDS" in result.output

    def test_empty_gff_raises_error(self, runner, ref_fasta, vcf_file,
                                     tmp_path):
        """Completely empty GFF (header only) → exit 1."""
        gff_empty = tmp_path / "empty.gff3"
        gff_empty.write_text("##gff-version 3\n")
        result = runner.invoke(main, [
            "run", "--vcf", vcf_file, "--gff", str(gff_empty),
            "--fasta", ref_fasta, "--outdir", str(tmp_path / "out"),
            "--min-freq", "0", "--min-depth", "0", "--min-qual", "0",
        ])
        assert result.exit_code == 1

    def test_run_parallel_raises_NoGenesFoundError(self, ref_fasta, vcf_file,
                                                    tmp_path):
        """Direct call to run_parallel raises NoGenesFoundError."""
        gff_no_cds = tmp_path / "no_cds.gff3"
        gff_no_cds.write_text(textwrap.dedent("""\
            ##gff-version 3
            chr1\ttest\tgene\t1\t90\t.\t+\t.\tID=gene1
            chr1\ttest\tmRNA\t1\t90\t.\t+\t.\tID=mRNA1;Parent=gene1
        """))
        with pytest.raises(NoGenesFoundError, match="No genes with CDS"):
            run_parallel(ref_fasta, str(gff_no_cds), vcf_file,
                         min_freq=0, min_depth=0, min_qual=0)


# ---------------------------------------------------------------------------
# 2. Zero shared contigs → ValueError fail-fast
# ---------------------------------------------------------------------------

class TestZeroSharedContigs:
    def test_completely_mismatched_contigs_raises(self, ref_fasta, gff3_file,
                                                   tmp_path):
        """VCF contigs (scaffold_1) vs GFF contigs (chr1) → ValueError."""
        vcf_wrong = tmp_path / "wrong_contigs.vcf"
        vcf_wrong.write_text(textwrap.dedent("""\
            ##fileformat=VCFv4.2
            ##contig=<ID=scaffold_1,length=350>
            #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
            scaffold_1\t10\t.\tA\tG\t50\t.\t.
        """))
        vcf_gz = _bgzip_and_index(vcf_wrong)

        with pytest.raises(ValueError, match="No contig names shared"):
            run_parallel(ref_fasta, gff3_file, vcf_gz,
                         min_freq=0, min_depth=0, min_qual=0)

    def test_cli_exits_with_error(self, runner, ref_fasta, gff3_file,
                                   tmp_path):
        """CLI catches ValueError and exits 1 with helpful message."""
        vcf_wrong = tmp_path / "wrong_contigs.vcf"
        vcf_wrong.write_text(textwrap.dedent("""\
            ##fileformat=VCFv4.2
            ##contig=<ID=scaffold_1,length=350>
            #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
            scaffold_1\t10\t.\tA\tG\t50\t.\t.
        """))
        vcf_gz = _bgzip_and_index(vcf_wrong)

        result = runner.invoke(main, [
            "run", "--vcf", vcf_gz, "--gff", gff3_file,
            "--fasta", ref_fasta, "--outdir", str(tmp_path / "out"),
            "--min-freq", "0", "--min-depth", "0", "--min-qual", "0",
        ])
        assert result.exit_code == 1
        assert "No contig names shared" in result.output


# ---------------------------------------------------------------------------
# 3. Partial contig mismatch → missing-contig genes dropped
# ---------------------------------------------------------------------------

class TestPartialContigMismatch:
    """Genes on contigs absent from VCF should be excluded from results."""

    @pytest.fixture
    def two_contig_data(self, tmp_path):
        """Create a dataset with genes on chr1 and chr2, VCF only has chr1."""
        # Reference: chr1 (350bp) + chr2 (90bp)
        ref_path = _write_fasta(tmp_path / "two_contig.fa", {
            "chr1": (
                "ATGGCTGATGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTG"
                "CTGCTGCTGCTGCTGCTTAAAAAAAAAAAAATGGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTG"
                "CTGCTGCTGCTGCTGCTACTTTTTTTTTTTTTTTTTTTTTTGCTGCTGCTGCTGATGCTGCTGCTGCTGC"
                "TGCTGCTTAAAAAAAAAAAATTAAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAG"
                "CAGCAGCAGCAGCAGCATCAGCAGCAGCCATAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
            ),
            "chr2": "ATGGCTGATGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTTAA",
        })

        # GFF: gene1 on chr1, gene_chr2 on chr2
        gff_path = tmp_path / "two_contig.gff3"
        gff_path.write_text(textwrap.dedent("""\
            ##gff-version 3
            chr1\ttest\tgene\t1\t90\t.\t+\t.\tID=gene1;Name=gene1
            chr1\ttest\tmRNA\t1\t90\t.\t+\t.\tID=mRNA1;Parent=gene1
            chr1\ttest\tCDS\t1\t90\t.\t+\t0\tID=cds1;Parent=mRNA1
            chr2\ttest\tgene\t1\t66\t.\t+\t.\tID=gene_chr2;Name=gene_chr2
            chr2\ttest\tmRNA\t1\t66\t.\t+\t.\tID=mRNA_chr2;Parent=gene_chr2
            chr2\ttest\tCDS\t1\t66\t.\t+\t0\tID=cds_chr2;Parent=mRNA_chr2
        """))

        # VCF: only chr1
        vcf_path = tmp_path / "chr1_only.vcf"
        vcf_path.write_text(textwrap.dedent("""\
            ##fileformat=VCFv4.2
            ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
            ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
            ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
            ##contig=<ID=chr1,length=350>
            #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
            chr1\t6\t.\tT\tC\t30\t.\t.\tGT:DP:AD\t0/1:100:80,20
        """))
        vcf_gz = _bgzip_and_index(vcf_path)

        return ref_path, str(gff_path), vcf_gz

    def test_missing_contig_genes_excluded(self, two_contig_data):
        """Genes on chr2 (absent from VCF) should not appear in results."""
        ref_path, gff_path, vcf_gz = two_contig_data
        results = run_parallel(ref_path, gff_path, vcf_gz,
                               min_freq=0, min_depth=0, min_qual=0)

        gene_ids = [r.gene_id for r in results]
        assert "gene1" in gene_ids
        assert "gene_chr2" not in gene_ids

    def test_missing_contig_genes_do_not_bias_summary(self, runner,
                                                        two_contig_data,
                                                        tmp_path):
        """Genome-wide piN/piS should not include zero-variant chr2 genes."""
        ref_path, gff_path, vcf_gz = two_contig_data
        outdir = tmp_path / "out"
        result = runner.invoke(main, [
            "run", "--vcf", vcf_gz, "--gff", gff_path,
            "--fasta", ref_path, "--outdir", str(outdir),
            "--min-freq", "0", "--min-depth", "0", "--min-qual", "0",
        ])
        assert result.exit_code == 0, result.output

        gene_df = pd.read_csv(outdir / "gene_results.tsv", sep="\t")
        assert len(gene_df) == 1  # only chr1 gene
        assert gene_df.iloc[0]["gene_id"] == "gene1"

        summary_df = pd.read_csv(outdir / "summary.tsv", sep="\t")
        assert summary_df.iloc[0]["total_genes"] == 1


# ---------------------------------------------------------------------------
# 4 & 5. Ambiguous reference codon (N bases) filtering
# ---------------------------------------------------------------------------

class TestAmbiguousReferenceCodon:
    """Reference genome with N bases should not crash; those codons are skipped."""

    @pytest.fixture
    def ref_with_n(self, tmp_path):
        """Reference with N bases in the middle of a CDS."""
        # 90bp: ATG GCT NNN GCT × 26 TAA
        # codon 3 (NNN) should be skipped
        seq = "ATGGCT" + "NNN" + "GCT" * 26 + "TAA"
        return _write_fasta(tmp_path / "ref_with_n.fa", {"chr1": seq})

    @pytest.fixture
    def gff_one_gene(self, tmp_path):
        gff = tmp_path / "one_gene.gff3"
        gff.write_text(textwrap.dedent("""\
            ##gff-version 3
            chr1\ttest\tgene\t1\t90\t.\t+\t.\tID=gene1
            chr1\ttest\tmRNA\t1\t90\t.\t+\t.\tID=mRNA1;Parent=gene1
            chr1\ttest\tCDS\t1\t90\t.\t+\t0\tID=cds1;Parent=mRNA1
        """))
        return str(gff)

    @pytest.fixture
    def empty_vcf(self, tmp_path):
        """VCF with chr1 contig but no variant records."""
        vcf = tmp_path / "empty.vcf"
        vcf.write_text(textwrap.dedent("""\
            ##fileformat=VCFv4.2
            ##contig=<ID=chr1,length=90>
            #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
        """))
        return _bgzip_and_index(vcf)

    def test_n_bases_skipped_without_crash(self, runner, ref_with_n,
                                            gff_one_gene, empty_vcf, tmp_path):
        """Pipeline completes successfully when reference has N bases."""
        result = runner.invoke(main, [
            "run", "--vcf", empty_vcf, "--gff", gff_one_gene,
            "--fasta", ref_with_n, "--outdir", str(tmp_path / "out"),
            "--min-freq", "0", "--min-depth", "0", "--min-qual", "0",
        ])
        assert result.exit_code == 0, result.output

        df = pd.read_csv(tmp_path / "out" / "gene_results.tsv", sep="\t")
        g = df.iloc[0]
        # 30 codons total, minus 1 stop, minus 1 NNN = 28 analyzed
        assert g["n_codons"] == 28
        assert g["N_sites"] > 0
        assert g["n_variants"] == 0

    def test_all_n_reference_returns_empty(self, tmp_path):
        """CDS with all-N reference → GeneResult with n_codons=0."""
        from pie.annotation import GeneModel
        from pie.reference import ReferenceGenome
        from pie.vcf import VariantReader

        seq = "N" * 90
        ref_path = _write_fasta(tmp_path / "all_n.fa", {"chr1": seq})

        gff = tmp_path / "gene.gff3"
        gff.write_text(textwrap.dedent("""\
            ##gff-version 3
            chr1\ttest\tgene\t1\t90\t.\t+\t.\tID=gene1
            chr1\ttest\tmRNA\t1\t90\t.\t+\t.\tID=mRNA1;Parent=gene1
            chr1\ttest\tCDS\t1\t90\t.\t+\t0\tID=cds1;Parent=mRNA1
        """))

        vcf = tmp_path / "empty.vcf"
        vcf.write_text(textwrap.dedent("""\
            ##fileformat=VCFv4.2
            ##contig=<ID=chr1,length=90>
            #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
        """))
        vcf_gz = _bgzip_and_index(vcf)

        from pie.annotation import parse_annotations
        genes = parse_annotations(str(gff))
        ref = ReferenceGenome(ref_path)
        vcf_reader = VariantReader(vcf_gz, min_freq=0, min_depth=0, min_qual=0)

        result = compute_gene_diversity(genes[0], ref, vcf_reader)
        ref.close()
        vcf_reader.close()

        assert result.n_codons == 0
        assert result.N_sites == 0.0
        assert result.S_sites == 0.0
        assert result.n_variants == 0


# ---------------------------------------------------------------------------
# 6. summary.tsv uses cds_snp_variants column name
# ---------------------------------------------------------------------------

class TestCdsSnpVariantsColumn:
    def test_summary_column_is_cds_snp_variants(self, runner, ref_fasta,
                                                  gff3_file, vcf_file,
                                                  tmp_path):
        """summary.tsv should have cds_snp_variants, not total_variants."""
        result = runner.invoke(main, [
            "run", "--vcf", vcf_file, "--gff", gff3_file,
            "--fasta", ref_fasta, "--outdir", str(tmp_path),
            "--min-freq", "0", "--min-depth", "0", "--min-qual", "0",
        ])
        assert result.exit_code == 0, result.output

        df = pd.read_csv(tmp_path / "summary.tsv", sep="\t")
        assert "cds_snp_variants" in df.columns
        assert "total_variants" not in df.columns
