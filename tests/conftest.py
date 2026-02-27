# tests/conftest.py
import pysam
import pytest
from pathlib import Path

DATA_DIR = Path(__file__).parent / "data"
REAL_DATA_DIR = Path(__file__).parents[1] / "data" / "Acyrthosiphon_pisum"


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

def _bgzip_and_index(vcf_path):
    """Bgzip and tabix-index a plain VCF file. Returns .vcf.gz path."""
    gz_path = str(vcf_path) + ".gz"
    pysam.tabix_compress(str(vcf_path), gz_path, force=True)
    pysam.tabix_index(gz_path, preset="vcf", force=True)
    return gz_path


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
    return _bgzip_and_index(vcf_path)


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
    return _bgzip_and_index(vcf_path)


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
    return _bgzip_and_index(vcf_path)


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
    return _bgzip_and_index(vcf_path)


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
