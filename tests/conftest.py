# tests/conftest.py
import pytest
from pathlib import Path

DATA_DIR = Path(__file__).parent / "data"


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
