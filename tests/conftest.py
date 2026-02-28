# tests/conftest.py
import pytest
from click.testing import CliRunner
from pathlib import Path

from tests.helpers import bgzip_and_index, write_fasta  # noqa: F401 (re-exported)

DATA_DIR = Path(__file__).parent / "data"
REAL_DATA_DIR = Path(__file__).parents[1] / "data" / "Acyrthosiphon_pisum"


@pytest.fixture
def runner():
    return CliRunner()


@pytest.fixture
def ref_fasta():
    return str(DATA_DIR / "ref.fa")


@pytest.fixture
def gff3_file():
    return str(DATA_DIR / "genes.gff3")


@pytest.fixture
def gtf_file():
    return str(DATA_DIR / "genes.gtf")


@pytest.fixture
def vcf_file():
    return str(DATA_DIR / "variants.vcf.gz")


@pytest.fixture
def plain_vcf_file():
    return str(DATA_DIR / "variants.vcf")


# --- Real dataset fixtures (Acyrthosiphon pisum, 400 genes) ---

@pytest.fixture
def multiallelic_vcf_file(tmp_path):
    """VCF with a multiallelic site (two decomposed records at pos 7)."""
    vcf_content = """\
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##contig=<ID=chr1,length=350>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
chr1\t6\t.\tT\tC\t30\t.\t.\tGT:DP:AD\t0/1:100:80,20
chr1\t7\t.\tG\tA\t50\t.\t.\tGT:DP:AD\t0/1:100:50,30
chr1\t7\t.\tG\tC\t50\t.\t.\tGT:DP:AD\t0/1:100:50,20
chr1\t195\t.\tA\tT\t45\t.\t.\tGT:DP:AD\t0/1:100:60,40
"""
    vcf_path = tmp_path / "multiallelic.vcf"
    vcf_path.write_text(vcf_content)
    return bgzip_and_index(vcf_path)


@pytest.fixture
def multiallelic_inline_vcf_file(tmp_path):
    """VCF with a non-decomposed multiallelic site (single line, ALT=A,C at pos 7)."""
    vcf_content = """\
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##contig=<ID=chr1,length=350>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
chr1\t6\t.\tT\tC\t30\t.\t.\tGT:DP:AD\t0/1:100:80,20
chr1\t7\t.\tG\tA,C\t50\t.\t.\tGT:DP:AD\t0/1:100:50,30,20
chr1\t195\t.\tA\tT\t45\t.\t.\tGT:DP:AD\t0/1:100:60,40
"""
    vcf_path = tmp_path / "multiallelic_inline.vcf"
    vcf_path.write_text(vcf_content)
    return bgzip_and_index(vcf_path)


@pytest.fixture
def vcf_no_ad_info_af(tmp_path):
    """VCF without FORMAT/AD — uses INFO/DP + INFO/AF fallback."""
    vcf_content = """\
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##contig=<ID=chr1,length=350>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
chr1\t6\t.\tT\tC\t30\t.\tDP=100;AF=0.20\tGT:DP\t0/1:100
chr1\t195\t.\tA\tT\t45\t.\tDP=100;AF=0.40\tGT:DP\t0/1:100
"""
    vcf_path = tmp_path / "no_ad_info_af.vcf"
    vcf_path.write_text(vcf_content)
    return bgzip_and_index(vcf_path)


@pytest.fixture
def vcf_no_ad_no_format_dp(tmp_path):
    """VCF without FORMAT/AD and without FORMAT/DP — uses INFO/DP only."""
    vcf_content = """\
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=chr1,length=350>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
chr1\t6\t.\tT\tC\t30\t.\tDP=80;AF=0.25\tGT\t0/1
chr1\t195\t.\tA\tT\t45\t.\tDP=120;AF=0.50\tGT\t0/1
"""
    vcf_path = tmp_path / "no_ad_no_format_dp.vcf"
    vcf_path.write_text(vcf_content)
    return bgzip_and_index(vcf_path)


@pytest.fixture
def vcf_multiallelic_no_ad(tmp_path):
    """Multiallelic VCF without FORMAT/AD — keep-mode fallback uses original freq."""
    vcf_content = """\
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##contig=<ID=chr1,length=350>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
chr1\t6\t.\tT\tC\t30\t.\tDP=100;AF=0.20\tGT:DP\t0/1:100
chr1\t7\t.\tG\tA\t50\t.\tDP=100;AF=0.30\tGT:DP\t0/1:100
chr1\t7\t.\tG\tC\t50\t.\tDP=100;AF=0.20\tGT:DP\t0/1:100
chr1\t195\t.\tA\tT\t45\t.\tDP=100;AF=0.40\tGT:DP\t0/1:100
"""
    vcf_path = tmp_path / "multiallelic_no_ad.vcf"
    vcf_path.write_text(vcf_content)
    return bgzip_and_index(vcf_path)


@pytest.fixture
def vcf_no_ad_no_af(tmp_path):
    """VCF without FORMAT/AD and without INFO/AF — freq should be 0.0."""
    vcf_content = """\
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##contig=<ID=chr1,length=350>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
chr1\t6\t.\tT\tC\t30\t.\tDP=100\tGT:DP\t0/1:100
"""
    vcf_path = tmp_path / "no_ad_no_af.vcf"
    vcf_path.write_text(vcf_content)
    return bgzip_and_index(vcf_path)


@pytest.fixture
def individual_vcf_file(tmp_path):
    """Multi-sample VCF with GT fields for individual-mode testing.

    4 diploid samples (S1-S4), same positions as test data:
      pos 6   T>C  S1:0/1  S2:0/0  S3:0/1  S4:./.  -> called=3, AN=6, AC=2, freq=1/3, call_rate=0.75
      pos 7   G>A  S1:0/0  S2:0/1  S3:0/1  S4:0/0  -> called=4, AN=8, AC=2, freq=1/4, call_rate=1.00
      pos 195 A>T  S1:0/1  S2:0/1  S3:1/1  S4:0/1  -> called=4, AN=8, AC=5, freq=5/8, call_rate=1.00
      pos 297 A>G  S1:0/1  S2:./.  S3:./.  S4:0/0  -> called=2, AN=4, AC=1, freq=1/4, call_rate=0.50
    """
    vcf_content = """\
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=chr1,length=350>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4
chr1\t6\t.\tT\tC\t30\t.\t.\tGT\t0/1\t0/0\t0/1\t./.
chr1\t7\t.\tG\tA\t50\t.\t.\tGT\t0/0\t0/1\t0/1\t0/0
chr1\t195\t.\tA\tT\t45\t.\t.\tGT\t0/1\t0/1\t1/1\t0/1
chr1\t297\t.\tA\tG\t15\t.\t.\tGT\t0/1\t./.\t./.\t0/0
"""
    vcf_path = tmp_path / "individual.vcf"
    vcf_path.write_text(vcf_content)
    return bgzip_and_index(vcf_path)


@pytest.fixture
def individual_multiallelic_vcf_file(tmp_path):
    """Multi-sample VCF with a multiallelic site for individual-mode testing.

      pos 6   T>C    S1:0/1  S2:0/0  S3:0/1  S4:0/0  -> AN=8, AC=2, freq=1/4
      pos 7   G>A,C  S1:0/1  S2:0/2  S3:1/2  S4:1/2  -> AN=8, AC_A=3, AC_C=3, freq_A=3/8, freq_C=3/8
      pos 195 A>T    S1:0/1  S2:0/0  S3:0/1  S4:0/1  -> AN=8, AC=3, freq=3/8
    """
    vcf_content = """\
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=chr1,length=350>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4
chr1\t6\t.\tT\tC\t30\t.\t.\tGT\t0/1\t0/0\t0/1\t0/0
chr1\t7\t.\tG\tA,C\t50\t.\t.\tGT\t0/1\t0/2\t1/2\t1/2
chr1\t195\t.\tA\tT\t45\t.\t.\tGT\t0/1\t0/0\t0/1\t0/1
"""
    vcf_path = tmp_path / "individual_multiallelic.vcf"
    vcf_path.write_text(vcf_content)
    return bgzip_and_index(vcf_path)


# --- Robustness fixtures (PR#9: Issues #2, #4, #7) ---

@pytest.fixture
def ref_with_n_fasta(tmp_path):
    """Reference FASTA with N bases in coding regions (Issue #2).

    Same layout as ref.fa but with N bases injected into gene1 codons:
      codon 1: ATG (start) - kept clean
      codon 2: NCT        - N at position 4 (0-based 3)
      codon 3: GAT        - clean
      codons 4-30: GCT*27 - clean
    Gene1 occupies pos 1-90 (1-based), gene2/gene3 same as standard.
    """
    # Start with the standard ref.fa sequence, inject N at position 3 (0-based)
    seq = list(
        "ATGGCTGATGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTG"
        "CTGCTGCTGCTGCTGCTTAAAAAAAAAAAAATGGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTG"
        "CTGCTGCTGCTGCTGCTACTTTTTTTTTTTTTTTTTTTTTTGCTGCTGCTGCTGATGCTGCTGCTGCTGC"
        "TGCTGCTTAAAAAAAAAAAATTAAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAG"
        "CAGCAGCAGCAGCAGCATCAGCAGCAGCCATAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
    )
    seq[3] = "N"  # codon 2 becomes NCT
    return write_fasta(tmp_path / "ref_with_n.fa", {"chr1": "".join(seq)})


@pytest.fixture
def ref_all_n_fasta(tmp_path):
    """Reference FASTA where gene1 region is entirely N (edge case for Issue #2).

    Gene1 (pos 1-90) is all N's, rest is normal.
    """
    seq = list(
        "ATGGCTGATGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTG"
        "CTGCTGCTGCTGCTGCTTAAAAAAAAAAAAATGGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTG"
        "CTGCTGCTGCTGCTGCTACTTTTTTTTTTTTTTTTTTTTTTGCTGCTGCTGCTGATGCTGCTGCTGCTGC"
        "TGCTGCTTAAAAAAAAAAAATTAAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAG"
        "CAGCAGCAGCAGCAGCATCAGCAGCAGCCATAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
    )
    for i in range(90):
        seq[i] = "N"
    return write_fasta(tmp_path / "ref_all_n.fa", {"chr1": "".join(seq)})


@pytest.fixture
def cdsonly_gff(tmp_path):
    """GFF3 with only CDS/exon features, no gene parent (Issue #4).

    parse_annotations requires 'gene' features; this GFF should yield 0 genes.
    """
    content = """\
##gff-version 3
##sequence-region chr1 1 350
chr1\ttest\texon\t1\t90\t.\t+\t.\tID=exon1
chr1\ttest\tCDS\t1\t90\t.\t+\t0\tID=cds1
chr1\ttest\texon\t101\t220\t.\t+\t.\tID=exon2
chr1\ttest\tCDS\t101\t220\t.\t+\t0\tID=cds2
"""
    gff_path = tmp_path / "cdsonly.gff3"
    gff_path.write_text(content)
    return str(gff_path)


@pytest.fixture
def mismatch_vcf_file(tmp_path):
    """VCF with bare contig names ('1') instead of 'chr1' (Issue #7).

    Same variants as standard test data but contig renamed from chr1 -> 1.
    """
    vcf_content = """\
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##contig=<ID=1,length=350>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
1\t6\t.\tT\tC\t30\t.\t.\tGT:DP:AD\t0/1:100:80,20
1\t7\t.\tG\tA\t50\t.\t.\tGT:DP:AD\t0/1:100:70,30
1\t195\t.\tA\tT\t45\t.\t.\tGT:DP:AD\t0/1:100:60,40
1\t297\t.\tA\tG\t15\t.\t.\tGT:DP:AD\t0/1:100:50,50
"""
    vcf_path = tmp_path / "mismatch.vcf"
    vcf_path.write_text(vcf_content)
    return bgzip_and_index(vcf_path)


@pytest.fixture
def real_ref_fasta():
    path = REAL_DATA_DIR / "Acyrthosiphon_pisum.fa"
    if not path.exists():
        pytest.skip("Real test data not extracted (run: tar xzf data/Acyrthosiphon_pisum.tar.gz -C data/)")
    return str(path)


@pytest.fixture
def real_gff():
    path = REAL_DATA_DIR / "Acyrthosiphon_pisum.gff"
    if not path.exists():
        pytest.skip("Real test data not extracted (run: tar xzf data/Acyrthosiphon_pisum.tar.gz -C data/)")
    return str(path)


@pytest.fixture
def real_vcf():
    path = REAL_DATA_DIR / "SRR27175631.filtered.snps.vcf.gz"
    if not path.exists():
        pytest.skip("Real test data not extracted (run: tar xzf data/Acyrthosiphon_pisum.tar.gz -C data/)")
    return str(path)
