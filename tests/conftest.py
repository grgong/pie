# tests/conftest.py
import subprocess

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
    with open(gz_path, "wb") as out:
        subprocess.run(["bgzip", "-c", str(vcf_path)], stdout=out, check=True)
    subprocess.run(["tabix", "-p", "vcf", gz_path], check=True)
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
