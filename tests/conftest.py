# tests/conftest.py
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
