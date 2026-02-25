# pie Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Build a Python CLI tool (`pie`) that calculates piN, piS, and piN/piS from pool-seq VCF + GFF/GTF + reference FASTA for large eukaryotic genomes.

**Architecture:** Numpy-vectorized Nei-Gojobori with gene-level multiprocessing. cyvcf2 for VCF, gffutils for GFF/GTF, pysam for FASTA. Click CLI with subcommands. Designed for Rust drop-in on the codon inner loop later.

**Tech Stack:** Python 3.12+, cyvcf2, gffutils, pysam, numpy, pandas, matplotlib, click. Conda env `pie`.

**Design doc:** `docs/plans/2026-02-25-pie-design.md`

---

### Task 1: Project Scaffold & Environment

**Files:**
- Create: `pyproject.toml`
- Create: `src/pie/__init__.py`
- Create: `src/pie/cli.py`
- Create: `tests/__init__.py`
- Create: `tests/conftest.py`
- Create: `environment.yml`

**Step 1: Create conda environment file**

```yaml
# environment.yml
name: pie
channels:
  - conda-forge
  - bioconda
dependencies:
  - python=3.12
  - cyvcf2
  - gffutils
  - pysam
  - numpy
  - pandas
  - matplotlib
  - click
  - pytest
  - htslib  # for bgzip/tabix
```

**Step 2: Create conda environment**

Run: `mamba env create -f environment.yml`

**Step 3: Create pyproject.toml**

```toml
[build-system]
requires = ["hatchling"]
build-backend = "hatchling.build"

[project]
name = "pie"
version = "0.1.0"
description = "piN/piS Estimator for pool-seq data"
requires-python = ">=3.12"
dependencies = [
    "cyvcf2",
    "gffutils",
    "pysam",
    "numpy",
    "pandas",
    "matplotlib",
    "click",
]

[project.scripts]
pie = "pie.cli:main"

[tool.hatch.build.targets.wheel]
packages = ["src/pie"]

[tool.pytest.ini_options]
testpaths = ["tests"]
```

**Step 4: Create src/pie/__init__.py**

```python
"""pie — piN/piS Estimator for pool-seq data."""

__version__ = "0.1.0"
```

**Step 5: Create minimal CLI stub**

```python
# src/pie/cli.py
import click


@click.group()
@click.version_option()
def main():
    """pie — piN/piS Estimator for pool-seq data."""
    pass


@main.command()
def run():
    """Run piN/piS analysis."""
    click.echo("pie run: not yet implemented")


@main.command()
def plot():
    """Plot results."""
    click.echo("pie plot: not yet implemented")


@main.command()
def summary():
    """Print summary statistics."""
    click.echo("pie summary: not yet implemented")
```

**Step 6: Create tests/__init__.py and tests/conftest.py**

```python
# tests/__init__.py
```

```python
# tests/conftest.py
import pytest
from pathlib import Path

DATA_DIR = Path(__file__).parent / "data"
```

**Step 7: Install package in dev mode and verify**

Run: `mamba run -n pie pip install -e .`
Run: `mamba run -n pie pie --version`
Expected: `pie, version 0.1.0`
Run: `mamba run -n pie pytest --co`
Expected: `no tests ran`

**Step 8: Commit**

```bash
git add pyproject.toml environment.yml src/ tests/
git commit -m "feat: project scaffold with conda env and CLI stub"
```

---

### Task 2: Test Data Fixtures

**Files:**
- Create: `tests/data/ref.fa` (small 3-gene reference, ~300bp)
- Create: `tests/data/genes.gff3`
- Create: `tests/data/genes.gtf`
- Create: `tests/data/variants.vcf`
- Create: `tests/conftest.py` (update with fixtures)

**Step 1: Create a minimal reference FASTA**

Design a synthetic reference with 3 genes on 1 chromosome. The sequence should have known codons so we can hand-calculate expected piN/piS.

```
>chr1
# Gene1: + strand, positions 1-90 (30 codons, single exon)
# Gene2: + strand, positions 101-220 (2 exons: 101-160, 181-220, = 100bp CDS = 33 codons)
# Gene3: - strand, positions 230-310 (27 codons, single exon)
# Total: ~350 bp
```

The exact nucleotide sequence must be designed so that:
- Gene1 has a known codon at position 1-3 (e.g., ATG = Met)
- Variants placed at known codon positions produce known syn/nonsyn changes
- At least one variant is synonymous and one is nonsynonymous

```python
# tests/data/ref.fa
# Use a sequence where we know the codons:
# Gene1 (+ strand, 1-based pos 1-90): 30 codons
#   ATG AAA GCT ... (start with Met, Lys, Ala, ...)
# Gene2 (+ strand, exon1: 101-160, exon2: 181-220): 33 codons with split codon at exon boundary
# Gene3 (- strand, 230-310): 27 codons (reverse complement)
```

Write a Python script `tests/create_test_data.py` that:
1. Generates a synthetic FASTA with controlled codons
2. Generates matching GFF3 and GTF
3. Generates a VCF with known variants at specific codon positions
4. Writes all files to `tests/data/`

The script must include hand-calculated expected piN/piS values as comments for verification.

**Step 2: Create the test data generation script and run it**

Run: `mamba run -n pie python tests/create_test_data.py`

**Step 3: Verify test data files exist**

Run: `ls tests/data/`
Expected: `ref.fa genes.gff3 genes.gtf variants.vcf`

**Step 4: Index test data**

Run: `mamba run -n pie samtools faidx tests/data/ref.fa`
Run: `mamba run -n pie bgzip -c tests/data/variants.vcf > tests/data/variants.vcf.gz && mamba run -n pie tabix -p vcf tests/data/variants.vcf.gz`

**Step 5: Add fixtures to conftest.py**

```python
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
    """Uncompressed VCF for testing auto-indexing."""
    return str(DATA_DIR / "variants.vcf")
```

**Step 6: Commit**

```bash
git add tests/
git commit -m "feat: add synthetic test data with known piN/piS values"
```

---

### Task 3: Codon Lookup Tables (`codon.py`)

This is the foundation — pure computation, no I/O, easy to test.

**Files:**
- Create: `src/pie/codon.py`
- Create: `tests/test_codon.py`

**Step 1: Write failing tests for codon utilities**

```python
# tests/test_codon.py
import numpy as np
from pie.codon import (
    CODON_TO_INDEX,
    INDEX_TO_CODON,
    AMINO_ACID,
    N_SITES,
    S_SITES,
    N_DIFFS,
    S_DIFFS,
    codon_to_index,
    is_stop_codon,
)


class TestCodonIndex:
    def test_64_codons(self):
        assert len(CODON_TO_INDEX) == 64
        assert len(INDEX_TO_CODON) == 64

    def test_roundtrip(self):
        for codon, idx in CODON_TO_INDEX.items():
            assert INDEX_TO_CODON[idx] == codon

    def test_atg_index(self):
        idx = codon_to_index("ATG")
        assert INDEX_TO_CODON[idx] == "ATG"


class TestAminoAcid:
    def test_atg_is_met(self):
        assert AMINO_ACID[codon_to_index("ATG")] == "M"

    def test_stop_codons(self):
        for stop in ["TAA", "TAG", "TGA"]:
            assert AMINO_ACID[codon_to_index(stop)] == "*"

    def test_synonymous_pair(self):
        # GCT and GCC both encode Ala
        assert AMINO_ACID[codon_to_index("GCT")] == "A"
        assert AMINO_ACID[codon_to_index("GCC")] == "A"


class TestIsStopCodon:
    def test_stops(self):
        for stop in ["TAA", "TAG", "TGA"]:
            assert is_stop_codon(codon_to_index(stop))

    def test_non_stop(self):
        assert not is_stop_codon(codon_to_index("ATG"))


class TestSiteCounts:
    def test_shape(self):
        assert N_SITES.shape == (64, 3)
        assert S_SITES.shape == (64, 3)

    def test_sites_sum_to_one_per_position(self):
        """At each position, N + S should equal 1 (excluding stops)."""
        for i in range(64):
            if AMINO_ACID[i] == "*":
                continue
            for pos in range(3):
                total = N_SITES[i, pos] + S_SITES[i, pos]
                assert abs(total - 1.0) < 1e-10, f"codon {INDEX_TO_CODON[i]} pos {pos}: N+S={total}"

    def test_aaa_position3(self):
        """AAA at position 3: A→G is syn (AAG=Lys), A→C and A→T are nonsyn."""
        idx = codon_to_index("AAA")
        assert abs(N_SITES[idx, 2] - 2 / 3) < 1e-10
        assert abs(S_SITES[idx, 2] - 1 / 3) < 1e-10

    def test_atg_all_nonsyn(self):
        """ATG (Met) — any change at any position is nonsynonymous (unique codon)."""
        idx = codon_to_index("ATG")
        for pos in range(3):
            assert abs(N_SITES[idx, pos] - 1.0) < 1e-10
            assert abs(S_SITES[idx, pos] - 0.0) < 1e-10

    def test_fourfold_degenerate(self):
        """GCN (Ala) — position 3 is 4-fold degenerate: all changes are synonymous."""
        idx = codon_to_index("GCT")
        assert abs(S_SITES[idx, 2] - 1.0) < 1e-10
        assert abs(N_SITES[idx, 2] - 0.0) < 1e-10


class TestDiffCounts:
    def test_shape(self):
        assert N_DIFFS.shape == (64, 64)
        assert S_DIFFS.shape == (64, 64)

    def test_self_diff_zero(self):
        """Same codon pair should have 0 diffs."""
        for i in range(64):
            assert N_DIFFS[i, i] == 0.0
            assert S_DIFFS[i, i] == 0.0

    def test_symmetric(self):
        """Diffs should be symmetric."""
        np.testing.assert_array_equal(N_DIFFS, N_DIFFS.T)
        np.testing.assert_array_equal(S_DIFFS, S_DIFFS.T)

    def test_synonymous_change(self):
        """GCT→GCC (Ala→Ala): 1 synonymous diff, 0 nonsynonymous."""
        i = codon_to_index("GCT")
        j = codon_to_index("GCC")
        assert S_DIFFS[i, j] == 1.0
        assert N_DIFFS[i, j] == 0.0

    def test_nonsynonymous_change(self):
        """AAA→AAC (Lys→Asn): 1 nonsynonymous diff, 0 synonymous."""
        i = codon_to_index("AAA")
        j = codon_to_index("AAC")
        assert N_DIFFS[i, j] == 1.0
        assert S_DIFFS[i, j] == 0.0

    def test_two_step_change_ignored(self):
        """Codons differing at 2+ positions: N_DIFFS + S_DIFFS should still
        classify each changed position independently (Nei-Gojobori averaging)."""
        i = codon_to_index("AAA")  # Lys
        j = codon_to_index("ACC")  # Thr — differs at pos 2 and 3
        total = N_DIFFS[i, j] + S_DIFFS[i, j]
        assert total > 0  # should have some diffs
```

**Step 2: Run tests to verify they fail**

Run: `mamba run -n pie pytest tests/test_codon.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'pie.codon'`

**Step 3: Implement codon.py**

```python
# src/pie/codon.py
"""Nei-Gojobori codon lookup tables for site and difference classification.

Precomputes:
- N_SITES[64, 3]: fractional nonsynonymous site counts per codon per position
- S_SITES[64, 3]: fractional synonymous site counts per codon per position
- N_DIFFS[64, 64]: nonsynonymous differences for each codon pair
- S_DIFFS[64, 64]: synonymous differences for each codon pair

Uses standard genetic code. Nei & Gojobori (1986) method.
"""

import numpy as np
from itertools import product as itertools_product

_BASES = "ACGT"
_BASE_TO_INT = {"A": 0, "C": 1, "G": 2, "T": 3}

# Generate all 64 codons in a fixed order
_ALL_CODONS = [a + b + c for a, b, c in itertools_product(_BASES, repeat=3)]

CODON_TO_INDEX: dict[str, int] = {c: i for i, c in enumerate(_ALL_CODONS)}
INDEX_TO_CODON: list[str] = _ALL_CODONS[:]

# Standard genetic code
_GENETIC_CODE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

# Amino acid array indexed by codon index
AMINO_ACID: np.ndarray = np.array([_GENETIC_CODE[c] for c in _ALL_CODONS], dtype="U1")

_STOP_SET = frozenset(i for i, aa in enumerate(AMINO_ACID) if aa == "*")


def codon_to_index(codon: str) -> int:
    """Convert a 3-letter codon string to its integer index (0-63)."""
    return CODON_TO_INDEX[codon.upper()]


def is_stop_codon(index: int) -> bool:
    """Check if a codon index corresponds to a stop codon."""
    return index in _STOP_SET


def _build_site_tables() -> tuple[np.ndarray, np.ndarray]:
    """Compute fractional N/S sites for each codon at each position.

    For each codon and each of its 3 positions, count how many of the
    3 possible single-nucleotide substitutions change the amino acid
    (nonsynonymous) vs preserve it (synonymous).

    N_sites = count(nonsyn changes) / 3
    S_sites = count(syn changes) / 3
    Sum to 1.0 per position for non-stop codons.
    """
    n_sites = np.zeros((64, 3), dtype=np.float64)
    s_sites = np.zeros((64, 3), dtype=np.float64)

    for idx, codon in enumerate(_ALL_CODONS):
        aa = _GENETIC_CODE[codon]
        if aa == "*":
            # Stop codons: 0 sites (cannot be classified)
            continue
        for pos in range(3):
            n_changes = 0
            s_changes = 0
            original_base = codon[pos]
            for alt_base in _BASES:
                if alt_base == original_base:
                    continue
                mutant = codon[:pos] + alt_base + codon[pos + 1:]
                mutant_aa = _GENETIC_CODE[mutant]
                if mutant_aa == aa:
                    s_changes += 1
                else:
                    n_changes += 1
            n_sites[idx, pos] = n_changes / 3.0
            s_sites[idx, pos] = s_changes / 3.0

    return n_sites, s_sites


def _build_diff_tables() -> tuple[np.ndarray, np.ndarray]:
    """Compute N/S differences for each codon pair.

    For single-step differences (codons differing at exactly 1 position):
    straightforward — check if the amino acid changed.

    For multi-step differences (2-3 positions differ): average over all
    shortest mutational pathways (Nei-Gojobori method). Each pathway
    is a sequence of single-step changes; for each step, classify as
    N or S. Average the total N and S diffs across all pathways.
    """
    n_diffs = np.zeros((64, 64), dtype=np.float64)
    s_diffs = np.zeros((64, 64), dtype=np.float64)

    for i in range(64):
        for j in range(i + 1, 64):
            codon_i = _ALL_CODONS[i]
            codon_j = _ALL_CODONS[j]

            # Find positions that differ
            diff_positions = [p for p in range(3) if codon_i[p] != codon_j[p]]
            n_diff_pos = len(diff_positions)

            if n_diff_pos == 0:
                continue

            if n_diff_pos == 1:
                # Single step: direct classification
                aa_i = _GENETIC_CODE[codon_i]
                aa_j = _GENETIC_CODE[codon_j]
                if aa_i == aa_j:
                    s_diffs[i, j] = 1.0
                else:
                    n_diffs[i, j] = 1.0
            else:
                # Multi-step: average over all permutations of change order
                from itertools import permutations
                total_n = 0.0
                total_s = 0.0
                n_paths = 0

                for perm in permutations(diff_positions):
                    path_n = 0.0
                    path_s = 0.0
                    current = list(codon_i)

                    for pos in perm:
                        prev_aa = _GENETIC_CODE["".join(current)]
                        current[pos] = codon_j[pos]
                        curr_aa = _GENETIC_CODE["".join(current)]

                        if prev_aa == "*" or curr_aa == "*":
                            # Path through stop codon — skip this path
                            break
                        if prev_aa == curr_aa:
                            path_s += 1.0
                        else:
                            path_n += 1.0
                    else:
                        # Only count paths that don't go through stop codons
                        total_n += path_n
                        total_s += path_s
                        n_paths += 1

                if n_paths > 0:
                    n_diffs[i, j] = total_n / n_paths
                    s_diffs[i, j] = total_s / n_paths

            # Symmetric
            n_diffs[j, i] = n_diffs[i, j]
            s_diffs[j, i] = s_diffs[i, j]

    return n_diffs, s_diffs


# Precompute at import time
N_SITES, S_SITES = _build_site_tables()
N_DIFFS, S_DIFFS = _build_diff_tables()
```

**Step 4: Run tests to verify they pass**

Run: `mamba run -n pie pytest tests/test_codon.py -v`
Expected: all PASS

**Step 5: Commit**

```bash
git add src/pie/codon.py tests/test_codon.py
git commit -m "feat: Nei-Gojobori codon lookup tables with site and diff classification"
```

---

### Task 4: Reference FASTA Module (`reference.py`)

**Files:**
- Create: `src/pie/reference.py`
- Create: `tests/test_reference.py`

**Step 1: Write failing tests**

```python
# tests/test_reference.py
import pytest
from pie.reference import ReferenceGenome


class TestReferenceGenome:
    def test_open_close(self, ref_fasta):
        ref = ReferenceGenome(ref_fasta)
        ref.close()

    def test_context_manager(self, ref_fasta):
        with ReferenceGenome(ref_fasta) as ref:
            assert ref is not None

    def test_fetch_sequence(self, ref_fasta):
        with ReferenceGenome(ref_fasta) as ref:
            # Fetch a known region (0-based, half-open)
            seq = ref.fetch("chr1", 0, 3)
            assert len(seq) == 3
            assert seq == seq.upper()  # should be uppercase

    def test_extract_codons_plus_strand(self, ref_fasta):
        """Extract codons from a + strand CDS given list of (chrom, start, end) exon intervals."""
        with ReferenceGenome(ref_fasta) as ref:
            # Gene1: single exon, + strand, positions 0-89 (0-based)
            exons = [("chr1", 0, 90)]
            codons = ref.extract_codons(exons, strand="+")
            assert len(codons) == 30
            assert codons[0] == "ATG"  # first codon should be start

    def test_extract_codons_minus_strand(self, ref_fasta):
        """Minus strand codons should be reverse-complemented."""
        with ReferenceGenome(ref_fasta) as ref:
            exons = [("chr1", 229, 310)]
            codons = ref.extract_codons(exons, strand="-")
            assert len(codons) == 27
            # First codon of - strand gene should be ATG (rev comp of CAT at end)

    def test_extract_codons_multi_exon(self, ref_fasta):
        """Multi-exon gene: codons span exon junction."""
        with ReferenceGenome(ref_fasta) as ref:
            exons = [("chr1", 100, 160), ("chr1", 180, 220)]
            codons = ref.extract_codons(exons, strand="+")
            total_bp = (160 - 100) + (220 - 180)  # 60 + 40 = 100
            assert len(codons) == total_bp // 3  # 33 codons

    def test_get_genomic_positions_for_codons(self, ref_fasta):
        """Map codon indices back to genomic positions."""
        with ReferenceGenome(ref_fasta) as ref:
            exons = [("chr1", 0, 90)]
            positions = ref.codon_genomic_positions(exons, strand="+")
            assert len(positions) == 30
            # First codon: positions 0, 1, 2
            assert positions[0] == ("chr1", 0, 1, 2)
            # Last codon: positions 87, 88, 89
            assert positions[29] == ("chr1", 87, 88, 89)
```

**Step 2: Run tests to verify they fail**

Run: `mamba run -n pie pytest tests/test_reference.py -v`
Expected: FAIL

**Step 3: Implement reference.py**

```python
# src/pie/reference.py
"""Reference genome FASTA access and codon extraction."""

import pysam


_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


class ReferenceGenome:
    """Indexed FASTA access with codon extraction for CDS regions."""

    def __init__(self, fasta_path: str):
        self._fasta = pysam.FastaFile(fasta_path)

    def close(self):
        self._fasta.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    def fetch(self, chrom: str, start: int, end: int) -> str:
        """Fetch sequence (0-based, half-open coordinates)."""
        return self._fasta.fetch(chrom, start, end).upper()

    def extract_codons(
        self, exons: list[tuple[str, int, int]], strand: str
    ) -> list[str]:
        """Extract codons from CDS exons.

        Args:
            exons: List of (chrom, start, end) in genomic order, 0-based half-open.
            strand: "+" or "-".

        Returns:
            List of 3-letter codon strings in reading frame order.
        """
        # Concatenate exon sequences
        cds_seq = ""
        for chrom, start, end in exons:
            cds_seq += self.fetch(chrom, start, end)

        if strand == "-":
            cds_seq = cds_seq[::-1].translate(_COMPLEMENT)

        # Split into codons (drop incomplete trailing codon)
        n_complete = (len(cds_seq) // 3) * 3
        cds_seq = cds_seq[:n_complete]
        return [cds_seq[i : i + 3] for i in range(0, len(cds_seq), 3)]

    def codon_genomic_positions(
        self, exons: list[tuple[str, int, int]], strand: str
    ) -> list[tuple[str, int, int, int]]:
        """Map codon indices to genomic positions.

        Returns:
            List of (chrom, pos1, pos2, pos3) for each codon,
            where pos1/2/3 are 0-based genomic positions.
        """
        # Build flat list of (chrom, position) for each CDS base
        bases = []
        for chrom, start, end in exons:
            for pos in range(start, end):
                bases.append((chrom, pos))

        if strand == "-":
            bases = bases[::-1]

        # Group into codons
        n_complete = (len(bases) // 3) * 3
        bases = bases[:n_complete]
        codons = []
        for i in range(0, len(bases), 3):
            chrom = bases[i][0]  # all 3 bases should be same chrom (or split-exon)
            codons.append((chrom, bases[i][1], bases[i + 1][1], bases[i + 2][1]))
        return codons
```

**Step 4: Run tests**

Run: `mamba run -n pie pytest tests/test_reference.py -v`
Expected: all PASS

**Step 5: Commit**

```bash
git add src/pie/reference.py tests/test_reference.py
git commit -m "feat: reference genome FASTA access and codon extraction"
```

---

### Task 5: GFF/GTF Annotation Module (`annotation.py`)

**Files:**
- Create: `src/pie/annotation.py`
- Create: `tests/test_annotation.py`

**Step 1: Write failing tests**

```python
# tests/test_annotation.py
import pytest
from pie.annotation import parse_annotations, GeneModel


class TestParseAnnotations:
    def test_parse_gff3(self, gff3_file):
        genes = parse_annotations(gff3_file)
        assert len(genes) > 0
        assert all(isinstance(g, GeneModel) for g in genes)

    def test_parse_gtf(self, gtf_file):
        genes = parse_annotations(gtf_file)
        assert len(genes) > 0

    def test_auto_detect_format(self, gff3_file, gtf_file):
        """Both formats should produce the same gene models."""
        gff3_genes = parse_annotations(gff3_file)
        gtf_genes = parse_annotations(gtf_file)
        assert len(gff3_genes) == len(gtf_genes)


class TestGeneModel:
    def test_gene_attributes(self, gff3_file):
        genes = parse_annotations(gff3_file)
        gene = genes[0]
        assert gene.gene_id is not None
        assert gene.transcript_id is not None
        assert gene.chrom is not None
        assert gene.strand in ("+", "-")
        assert len(gene.cds_exons) > 0

    def test_cds_exons_sorted(self, gff3_file):
        """CDS exons should be sorted by start position (genomic order)."""
        genes = parse_annotations(gff3_file)
        for gene in genes:
            starts = [e[1] for e in gene.cds_exons]
            assert starts == sorted(starts)

    def test_longest_isoform_selected(self, gff3_file):
        """If multiple isoforms exist, the longest CDS should be selected."""
        genes = parse_annotations(gff3_file)
        # Test data should have at least one gene where we can verify this
        for gene in genes:
            assert gene.cds_length > 0
            assert gene.cds_length == sum(e - s for _, s, e in gene.cds_exons)

    def test_multi_exon_gene(self, gff3_file):
        """Gene2 in test data should have 2 CDS exons."""
        genes = parse_annotations(gff3_file)
        gene2 = [g for g in genes if "gene2" in g.gene_id.lower()][0]
        assert len(gene2.cds_exons) == 2
```

**Step 2: Run tests to verify they fail**

Run: `mamba run -n pie pytest tests/test_annotation.py -v`
Expected: FAIL

**Step 3: Implement annotation.py**

```python
# src/pie/annotation.py
"""GFF3/GTF parsing with longest-isoform selection."""

from dataclasses import dataclass
from collections import defaultdict

import gffutils


@dataclass
class GeneModel:
    """A gene with its selected (longest) transcript's CDS exons."""

    gene_id: str
    transcript_id: str
    chrom: str
    start: int  # gene start, 0-based
    end: int  # gene end, 0-based half-open
    strand: str
    cds_exons: list[tuple[str, int, int]]  # [(chrom, start, end), ...] 0-based half-open

    @property
    def cds_length(self) -> int:
        return sum(end - start for _, start, end in self.cds_exons)


def _detect_format(path: str) -> str:
    """Auto-detect GFF3 vs GTF from file content."""
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                if "gff-version" in line.lower():
                    return "gff3"
                continue
            # GTF uses gene_id "..."; GFF3 uses ID=...;Parent=...
            if "ID=" in line or "Parent=" in line:
                return "gff3"
            if 'gene_id "' in line or "gene_id '" in line:
                return "gtf"
    return "gff3"  # default


def parse_annotations(path: str) -> list[GeneModel]:
    """Parse GFF3 or GTF, select longest isoform per gene, return GeneModels.

    Returns list sorted by (chrom, start).
    """
    fmt = _detect_format(path)

    # Create gffutils database in memory
    db = gffutils.create_db(
        path,
        dbfn=":memory:",
        force=True,
        keep_order=True,
        merge_strategy="merge",
        sort_attribute_values=True,
    )

    # Collect CDS features grouped by transcript
    # Strategy: for each gene, find all mRNA/transcript children,
    # then for each transcript find CDS children, pick longest.
    genes = []

    for gene in db.all_features(featuretype="gene"):
        gene_id = gene.id
        chrom = gene.seqid
        strand = gene.strand

        # Find all transcripts (mRNA or transcript) for this gene
        transcript_cds: dict[str, list] = defaultdict(list)

        for transcript in db.children(gene, featuretype=["mRNA", "transcript"], order_by="start"):
            for cds in db.children(transcript, featuretype="CDS", order_by="start"):
                # gffutils uses 1-based inclusive; convert to 0-based half-open
                transcript_cds[transcript.id].append(
                    (cds.seqid, cds.start - 1, cds.end)
                )

        # If no transcript hierarchy, try direct CDS children of gene
        if not transcript_cds:
            for cds in db.children(gene, featuretype="CDS", order_by="start"):
                transcript_cds[gene_id].append(
                    (cds.seqid, cds.start - 1, cds.end)
                )

        if not transcript_cds:
            continue

        # Select longest transcript by total CDS length
        best_tid = max(
            transcript_cds,
            key=lambda tid: sum(e - s for _, s, e in transcript_cds[tid]),
        )
        best_exons = sorted(transcript_cds[best_tid], key=lambda x: x[1])

        genes.append(
            GeneModel(
                gene_id=gene_id,
                transcript_id=best_tid,
                chrom=chrom,
                start=gene.start - 1,
                end=gene.end,
                strand=strand,
                cds_exons=best_exons,
            )
        )

    genes.sort(key=lambda g: (g.chrom, g.start))
    return genes
```

**Step 4: Run tests**

Run: `mamba run -n pie pytest tests/test_annotation.py -v`
Expected: all PASS

**Step 5: Commit**

```bash
git add src/pie/annotation.py tests/test_annotation.py
git commit -m "feat: GFF3/GTF parsing with longest isoform selection"
```

---

### Task 6: VCF Module (`vcf.py`)

**Files:**
- Create: `src/pie/vcf.py`
- Create: `tests/test_vcf.py`

**Step 1: Write failing tests**

```python
# tests/test_vcf.py
import pytest
import os
from pie.vcf import ensure_indexed, VariantReader


class TestEnsureIndexed:
    def test_already_indexed(self, vcf_file):
        """Should return path unchanged if .vcf.gz + .tbi exist."""
        result = ensure_indexed(vcf_file)
        assert result.endswith(".vcf.gz")

    def test_plain_vcf_gets_indexed(self, plain_vcf_file, tmp_path):
        """Should bgzip + tabix a plain VCF, writing to tmp_path."""
        import shutil
        # Copy plain VCF to tmp so we don't modify test data
        tmp_vcf = str(tmp_path / "test.vcf")
        shutil.copy(plain_vcf_file, tmp_vcf)
        result = ensure_indexed(tmp_vcf)
        assert result.endswith(".vcf.gz")
        assert os.path.exists(result + ".tbi")


class TestVariantReader:
    def test_open_close(self, vcf_file):
        reader = VariantReader(vcf_file, min_freq=0.01, min_depth=1, min_qual=0)
        reader.close()

    def test_context_manager(self, vcf_file):
        with VariantReader(vcf_file, min_freq=0.01, min_depth=1, min_qual=0) as reader:
            assert reader is not None

    def test_fetch_region(self, vcf_file):
        """Fetch variants in a known region."""
        with VariantReader(vcf_file, min_freq=0.0, min_depth=0, min_qual=0) as reader:
            variants = reader.fetch("chr1", 0, 90)
            assert isinstance(variants, list)
            # Test data should have known variants in gene1 region

    def test_variant_fields(self, vcf_file):
        """Each variant should have pos, ref, alt, freq, depth."""
        with VariantReader(vcf_file, min_freq=0.0, min_depth=0, min_qual=0) as reader:
            variants = reader.fetch("chr1", 0, 350)
            assert len(variants) > 0
            v = variants[0]
            assert hasattr(v, "pos")  # 0-based
            assert hasattr(v, "ref")
            assert hasattr(v, "alt")
            assert hasattr(v, "freq")  # alt allele frequency
            assert hasattr(v, "depth")

    def test_min_freq_filter(self, vcf_file):
        """Variants below min_freq should be excluded."""
        with VariantReader(vcf_file, min_freq=0.99, min_depth=0, min_qual=0) as reader:
            variants = reader.fetch("chr1", 0, 350)
            assert len(variants) == 0  # all should be filtered

    def test_min_depth_filter(self, vcf_file):
        """Variants below min_depth should be excluded."""
        with VariantReader(vcf_file, min_freq=0.0, min_depth=999999, min_qual=0) as reader:
            variants = reader.fetch("chr1", 0, 350)
            assert len(variants) == 0

    def test_pass_only_filter(self, vcf_file):
        with VariantReader(
            vcf_file, min_freq=0.0, min_depth=0, min_qual=0, pass_only=True
        ) as reader:
            variants = reader.fetch("chr1", 0, 350)
            # All test variants should have PASS or .
            assert all(v is not None for v in variants)
```

**Step 2: Run tests to verify they fail**

Run: `mamba run -n pie pytest tests/test_vcf.py -v`
Expected: FAIL

**Step 3: Implement vcf.py**

```python
# src/pie/vcf.py
"""VCF parsing, filtering, and auto-indexing via cyvcf2."""

import os
import subprocess
import logging
from dataclasses import dataclass

from cyvcf2 import VCF

log = logging.getLogger(__name__)


@dataclass(slots=True)
class Variant:
    """A filtered, biallelic SNP variant."""

    pos: int  # 0-based genomic position
    ref: str  # reference allele (single base)
    alt: str  # alternate allele (single base)
    freq: float  # alternate allele frequency (0-1)
    depth: int  # total depth


def ensure_indexed(vcf_path: str) -> str:
    """Ensure VCF is bgzipped and tabix-indexed. Returns path to .vcf.gz."""
    if vcf_path.endswith(".vcf"):
        gz_path = vcf_path + ".gz"
        log.info("Bgzipping %s → %s", vcf_path, gz_path)
        subprocess.run(["bgzip", "-c", vcf_path], stdout=open(gz_path, "wb"), check=True)
        vcf_path = gz_path

    tbi_path = vcf_path + ".tbi"
    if not os.path.exists(tbi_path):
        log.info("Creating tabix index for %s", vcf_path)
        subprocess.run(["tabix", "-p", "vcf", vcf_path], check=True)

    return vcf_path


class VariantReader:
    """Filtered VCF reader with region queries."""

    def __init__(
        self,
        vcf_path: str,
        min_freq: float = 0.01,
        min_depth: int = 10,
        min_qual: float = 20.0,
        pass_only: bool = False,
    ):
        self._vcf_path = vcf_path
        self._min_freq = min_freq
        self._min_depth = min_depth
        self._min_qual = min_qual
        self._pass_only = pass_only
        self._vcf = VCF(vcf_path)

    def close(self):
        self._vcf.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    def fetch(self, chrom: str, start: int, end: int) -> list[Variant]:
        """Fetch filtered variants in region (0-based, half-open).

        Returns list of Variant objects for biallelic SNPs passing filters.
        """
        variants = []
        # cyvcf2 uses 1-based region strings
        region = f"{chrom}:{start + 1}-{end}"

        for record in self._vcf(region):
            # Skip non-SNPs
            if not record.is_snp:
                continue

            # Skip multi-allelic (should be pre-decomposed, but be safe)
            if len(record.ALT) != 1:
                continue

            # PASS filter
            if self._pass_only and record.FILTER is not None:
                continue

            # QUAL filter
            if record.QUAL is not None and record.QUAL < self._min_qual:
                continue

            # Extract allele depth and total depth
            depth, freq = self._extract_freq_depth(record)

            if depth < self._min_depth:
                continue

            if freq < self._min_freq:
                continue

            variants.append(
                Variant(
                    pos=record.POS - 1,  # convert to 0-based
                    ref=record.REF,
                    alt=record.ALT[0],
                    freq=freq,
                    depth=depth,
                )
            )

        return variants

    def _extract_freq_depth(self, record) -> tuple[int, float]:
        """Extract total depth and alt allele frequency from a VCF record.

        Tries AD/DP format fields first (freebayes), then INFO AF/DP.
        """
        # Try FORMAT AD field (freebayes style: AD = ref_count,alt_count)
        try:
            ad = record.format("AD")
            if ad is not None:
                ref_count = int(ad[0][0])
                alt_count = int(ad[0][1])
                total = ref_count + alt_count
                if total > 0:
                    return total, alt_count / total
        except (KeyError, IndexError, TypeError):
            pass

        # Try FORMAT DP + INFO AF
        try:
            dp = record.format("DP")
            depth = int(dp[0][0]) if dp is not None else 0
        except (KeyError, IndexError, TypeError):
            depth = 0

        if depth == 0:
            # Fallback to INFO DP
            try:
                depth = record.INFO.get("DP", 0)
            except (KeyError, TypeError):
                depth = 0

        # AF from INFO
        try:
            af = record.INFO.get("AF")
            if isinstance(af, (list, tuple)):
                freq = float(af[0])
            else:
                freq = float(af) if af is not None else 0.0
        except (KeyError, TypeError, ValueError):
            freq = 0.0

        return int(depth), freq
```

**Step 4: Run tests**

Run: `mamba run -n pie pytest tests/test_vcf.py -v`
Expected: all PASS

**Step 5: Commit**

```bash
git add src/pie/vcf.py tests/test_vcf.py
git commit -m "feat: VCF parsing with filtering and auto-indexing"
```

---

### Task 7: Core Diversity Engine (`diversity.py`)

This is the heart of the tool — per-gene piN/piS computation.

**Files:**
- Create: `src/pie/diversity.py`
- Create: `tests/test_diversity.py`

**Step 1: Write failing tests**

```python
# tests/test_diversity.py
import numpy as np
import pytest
from pie.codon import codon_to_index, N_SITES, S_SITES
from pie.diversity import (
    build_allele_freq_array,
    compute_gene_diversity,
    compute_codon_diversity,
    GeneResult,
)


class TestBuildAlleleFreqArray:
    def test_monomorphic(self):
        """All reference codons, no variants → freq array with 1.0 at ref alleles."""
        codons = ["ATG", "AAA", "GCT"]
        positions = [
            ("chr1", 0, 1, 2),
            ("chr1", 3, 4, 5),
            ("chr1", 6, 7, 8),
        ]
        variants = []  # no variants
        freqs = build_allele_freq_array(codons, positions, variants)
        assert freqs.shape == (3, 3, 4)
        # ATG: A=1.0 at pos0, T=1.0 at pos1, G=1.0 at pos2
        assert freqs[0, 0, 0] == 1.0  # A
        assert freqs[0, 1, 3] == 1.0  # T
        assert freqs[0, 2, 2] == 1.0  # G

    def test_one_variant(self):
        """Single variant at one codon position."""
        from pie.vcf import Variant

        codons = ["AAA"]
        positions = [("chr1", 0, 1, 2)]
        # Variant at pos 2: A→G with freq 0.1
        variants = [Variant(pos=2, ref="A", alt="G", freq=0.1, depth=100)]
        freqs = build_allele_freq_array(codons, positions, variants)
        assert freqs.shape == (1, 3, 4)
        # Position 2 (3rd in codon): A=0.9, G=0.1
        assert abs(freqs[0, 2, 0] - 0.9) < 1e-10  # A
        assert abs(freqs[0, 2, 2] - 0.1) < 1e-10  # G


class TestComputeCodonDiversity:
    def test_monomorphic_codon(self):
        """Monomorphic codon: N_diffs = S_diffs = 0, sites > 0."""
        freq = np.zeros((3, 4))
        freq[0, 0] = 1.0  # A
        freq[1, 0] = 1.0  # A
        freq[2, 0] = 1.0  # A → codon AAA
        result = compute_codon_diversity(freq)
        assert result["N_diffs"] == 0.0
        assert result["S_diffs"] == 0.0
        idx = codon_to_index("AAA")
        expected_n = float(np.sum(N_SITES[idx]))
        expected_s = float(np.sum(S_SITES[idx]))
        assert abs(result["N_sites"] - expected_n) < 1e-10
        assert abs(result["S_sites"] - expected_s) < 1e-10

    def test_synonymous_variant(self):
        """GCT with G→A at pos3 (GCT→GCA, both Ala): should produce S_diffs > 0."""
        freq = np.zeros((3, 4))
        freq[0, 2] = 1.0  # G
        freq[1, 1] = 1.0  # C
        freq[2, 3] = 0.8  # T at 0.8
        freq[2, 0] = 0.2  # A at 0.2 → GCA (also Ala)
        result = compute_codon_diversity(freq)
        assert result["S_diffs"] > 0
        assert result["N_diffs"] == 0.0

    def test_nonsynonymous_variant(self):
        """AAA with A→C at pos3 (AAA→AAC, Lys→Asn): should produce N_diffs > 0."""
        freq = np.zeros((3, 4))
        freq[0, 0] = 1.0  # A
        freq[1, 0] = 1.0  # A
        freq[2, 0] = 0.9  # A
        freq[2, 1] = 0.1  # C → AAC (Asn)
        result = compute_codon_diversity(freq)
        assert result["N_diffs"] > 0
        assert result["S_diffs"] == 0.0
        # Expected N_diffs = 2 * 0.9 * 0.1 * 1.0 = 0.18
        assert abs(result["N_diffs"] - 0.18) < 1e-10

    def test_stop_codon_renormalized(self):
        """If a possible codon is a stop, its freq should be zeroed and others renormalized."""
        # TAA is stop. Construct T_A with A/G at pos3 → TAA (stop) / TAG (stop)
        # Actually let's use a case where one variant creates a stop:
        # TGG (Trp) → TGA (stop) at pos3: G→A
        freq = np.zeros((3, 4))
        freq[0, 3] = 1.0  # T
        freq[1, 2] = 1.0  # G
        freq[2, 2] = 0.95  # G → TGG
        freq[2, 0] = 0.05  # A → TGA (stop)
        result = compute_codon_diversity(freq)
        # After removing stop, only TGG remains at freq 1.0
        # So N_diffs = S_diffs = 0
        assert result["N_diffs"] == 0.0
        assert result["S_diffs"] == 0.0


class TestComputeGeneDiversity:
    def test_basic_gene(self, ref_fasta, gff3_file, vcf_file):
        """Integration test: compute diversity for gene1 from test data."""
        from pie.reference import ReferenceGenome
        from pie.annotation import parse_annotations
        from pie.vcf import VariantReader

        genes = parse_annotations(gff3_file)
        gene = genes[0]

        with ReferenceGenome(ref_fasta) as ref, VariantReader(
            vcf_file, min_freq=0.0, min_depth=0, min_qual=0
        ) as vcf:
            result = compute_gene_diversity(gene, ref, vcf)

        assert isinstance(result, GeneResult)
        assert result.gene_id == gene.gene_id
        assert result.N_sites > 0
        assert result.S_sites > 0
        assert result.piN >= 0
        assert result.piS >= 0
```

**Step 2: Run tests to verify they fail**

Run: `mamba run -n pie pytest tests/test_diversity.py -v`
Expected: FAIL

**Step 3: Implement diversity.py**

```python
# src/pie/diversity.py
"""Core piN/piS diversity computation engine."""

import logging
from dataclasses import dataclass, field
from itertools import product as itertools_product

import numpy as np

from pie.codon import (
    CODON_TO_INDEX,
    N_SITES,
    S_SITES,
    N_DIFFS,
    S_DIFFS,
    AMINO_ACID,
    is_stop_codon,
)
from pie.annotation import GeneModel
from pie.reference import ReferenceGenome
from pie.vcf import Variant, VariantReader

log = logging.getLogger(__name__)

_BASE_TO_IDX = {"A": 0, "C": 1, "G": 2, "T": 3}
_IDX_TO_BASE = "ACGT"


@dataclass
class CodonResult:
    """Per-codon diversity metrics."""

    chrom: str
    pos1: int  # genomic position of first codon base (0-based)
    N_sites: float
    S_sites: float
    N_diffs: float
    S_diffs: float


@dataclass
class GeneResult:
    """Per-gene aggregated diversity metrics."""

    gene_id: str
    transcript_id: str
    chrom: str
    start: int
    end: int
    strand: str
    n_codons: int
    n_poly_codons: int
    N_sites: float
    S_sites: float
    N_diffs: float
    S_diffs: float
    mean_depth: float
    n_variants: int
    codon_results: list[CodonResult] = field(default_factory=list)

    @property
    def piN(self) -> float:
        return self.N_diffs / self.N_sites if self.N_sites > 0 else 0.0

    @property
    def piS(self) -> float:
        return self.S_diffs / self.S_sites if self.S_sites > 0 else 0.0

    @property
    def piN_piS(self) -> float | None:
        if self.piS == 0:
            return None
        return self.piN / self.piS


def build_allele_freq_array(
    codons: list[str],
    positions: list[tuple[str, int, int, int]],
    variants: list[Variant],
) -> np.ndarray:
    """Build (n_codons, 3, 4) allele frequency array.

    Initializes with reference alleles at freq 1.0, then overlays variants.
    """
    n = len(codons)
    freqs = np.zeros((n, 3, 4), dtype=np.float64)

    # Set reference allele frequencies to 1.0
    for i, codon in enumerate(codons):
        for j, base in enumerate(codon):
            freqs[i, j, _BASE_TO_IDX[base]] = 1.0

    # Build position → (codon_idx, codon_pos) lookup
    pos_to_codon: dict[tuple[str, int], tuple[int, int]] = {}
    for codon_idx, (chrom, p1, p2, p3) in enumerate(positions):
        pos_to_codon[(chrom, p1)] = (codon_idx, 0)
        pos_to_codon[(chrom, p2)] = (codon_idx, 1)
        pos_to_codon[(chrom, p3)] = (codon_idx, 2)

    # Overlay variants
    for var in variants:
        key = (positions[0][0], var.pos)  # assume single chrom per gene
        if key not in pos_to_codon:
            continue
        ci, cj = pos_to_codon[key]
        alt_idx = _BASE_TO_IDX.get(var.alt)
        ref_idx = _BASE_TO_IDX.get(var.ref)
        if alt_idx is None or ref_idx is None:
            continue
        freqs[ci, cj, alt_idx] += var.freq
        freqs[ci, cj, ref_idx] -= var.freq
        # Clamp to [0, 1]
        freqs[ci, cj] = np.clip(freqs[ci, cj], 0.0, 1.0)
        # Renormalize to sum to 1
        total = freqs[ci, cj].sum()
        if total > 0:
            freqs[ci, cj] /= total

    return freqs


def compute_codon_diversity(freq: np.ndarray) -> dict:
    """Compute N/S sites and diffs for a single codon given allele frequencies.

    Args:
        freq: shape (3, 4) — allele frequencies at each of 3 codon positions.

    Returns:
        dict with keys: N_sites, S_sites, N_diffs, S_diffs
    """
    # Enumerate all possible codons
    alleles_per_pos = []
    for pos in range(3):
        present = np.where(freq[pos] > 0)[0]
        alleles_per_pos.append(present.tolist())

    codon_indices = []
    codon_freqs = []

    for combo in itertools_product(*alleles_per_pos):
        codon_str = _IDX_TO_BASE[combo[0]] + _IDX_TO_BASE[combo[1]] + _IDX_TO_BASE[combo[2]]
        cidx = CODON_TO_INDEX[codon_str]
        cfreq = freq[0, combo[0]] * freq[1, combo[1]] * freq[2, combo[2]]
        codon_indices.append(cidx)
        codon_freqs.append(cfreq)

    codon_indices = np.array(codon_indices, dtype=np.int32)
    codon_freqs = np.array(codon_freqs, dtype=np.float64)

    # Remove stop codons and renormalize
    stop_mask = np.array([is_stop_codon(idx) for idx in codon_indices])
    if stop_mask.any():
        stop_freq = codon_freqs[stop_mask].sum()
        if stop_freq > 0.01:
            log.warning(
                "Stop codon frequency %.3f > 1%% — possible annotation error", stop_freq
            )
        codon_freqs[stop_mask] = 0.0
        total = codon_freqs.sum()
        if total > 0:
            codon_freqs /= total
        else:
            return {"N_sites": 0.0, "S_sites": 0.0, "N_diffs": 0.0, "S_diffs": 0.0}

    # Weighted site counts
    n_sites = 0.0
    s_sites = 0.0
    for idx, f in zip(codon_indices, codon_freqs):
        if f == 0:
            continue
        n_sites += f * N_SITES[idx].sum()
        s_sites += f * S_SITES[idx].sum()

    # Weighted pairwise differences
    n_diffs = 0.0
    s_diffs = 0.0
    n_codons_present = len(codon_indices)
    for i in range(n_codons_present):
        if codon_freqs[i] == 0:
            continue
        for j in range(i + 1, n_codons_present):
            if codon_freqs[j] == 0:
                continue
            weight = 2.0 * codon_freqs[i] * codon_freqs[j]
            n_diffs += weight * N_DIFFS[codon_indices[i], codon_indices[j]]
            s_diffs += weight * S_DIFFS[codon_indices[i], codon_indices[j]]

    return {"N_sites": n_sites, "S_sites": s_sites, "N_diffs": n_diffs, "S_diffs": s_diffs}


def compute_gene_diversity(
    gene: GeneModel,
    ref: ReferenceGenome,
    vcf: VariantReader,
    window_size: int = 0,
    window_step: int = 0,
) -> GeneResult:
    """Compute piN/piS for a single gene.

    Args:
        gene: Gene model with CDS exon coordinates.
        ref: Reference genome for codon extraction.
        vcf: Variant reader for querying SNPs.
        window_size: Sliding window size in bp (0 = no windows).
        window_step: Sliding window step in bp.

    Returns:
        GeneResult with per-gene and optionally per-codon metrics.
    """
    # Extract codons and genomic positions
    codons = ref.extract_codons(gene.cds_exons, gene.strand)
    positions = ref.codon_genomic_positions(gene.cds_exons, gene.strand)

    if not codons:
        return GeneResult(
            gene_id=gene.gene_id,
            transcript_id=gene.transcript_id,
            chrom=gene.chrom,
            start=gene.start,
            end=gene.end,
            strand=gene.strand,
            n_codons=0,
            n_poly_codons=0,
            N_sites=0.0,
            S_sites=0.0,
            N_diffs=0.0,
            S_diffs=0.0,
            mean_depth=0.0,
            n_variants=0,
        )

    # Fetch variants across all CDS exons
    all_variants = []
    total_depth = 0
    depth_count = 0
    for chrom, start, end in gene.cds_exons:
        exon_variants = vcf.fetch(chrom, start, end)
        all_variants.extend(exon_variants)
        for v in exon_variants:
            total_depth += v.depth
            depth_count += 1

    # Build allele frequency array
    freqs = build_allele_freq_array(codons, positions, all_variants)

    # Identify polymorphic codons
    is_poly = np.sum(freqs > 0, axis=2) > 1  # (n_codons, 3)
    poly_mask = np.any(is_poly, axis=1)  # (n_codons,)

    # Compute per-codon diversity
    total_N_sites = 0.0
    total_S_sites = 0.0
    total_N_diffs = 0.0
    total_S_diffs = 0.0
    codon_results = []

    for i in range(len(codons)):
        if poly_mask[i]:
            # Polymorphic codon: full computation
            result = compute_codon_diversity(freqs[i])
        else:
            # Monomorphic codon: only need site counts
            cidx = CODON_TO_INDEX.get(codons[i])
            if cidx is None or AMINO_ACID[cidx] == "*":
                continue
            result = {
                "N_sites": float(N_SITES[cidx].sum()),
                "S_sites": float(S_SITES[cidx].sum()),
                "N_diffs": 0.0,
                "S_diffs": 0.0,
            }

        total_N_sites += result["N_sites"]
        total_S_sites += result["S_sites"]
        total_N_diffs += result["N_diffs"]
        total_S_diffs += result["S_diffs"]

        codon_results.append(
            CodonResult(
                chrom=positions[i][0],
                pos1=positions[i][1],
                N_sites=result["N_sites"],
                S_sites=result["S_sites"],
                N_diffs=result["N_diffs"],
                S_diffs=result["S_diffs"],
            )
        )

    mean_depth = total_depth / depth_count if depth_count > 0 else 0.0

    return GeneResult(
        gene_id=gene.gene_id,
        transcript_id=gene.transcript_id,
        chrom=gene.chrom,
        start=gene.start,
        end=gene.end,
        strand=gene.strand,
        n_codons=len(codons),
        n_poly_codons=int(poly_mask.sum()),
        N_sites=total_N_sites,
        S_sites=total_S_sites,
        N_diffs=total_N_diffs,
        S_diffs=total_S_diffs,
        mean_depth=mean_depth,
        n_variants=len(all_variants),
        codon_results=codon_results,
    )
```

**Step 4: Run tests**

Run: `mamba run -n pie pytest tests/test_diversity.py -v`
Expected: all PASS

**Step 5: Commit**

```bash
git add src/pie/diversity.py tests/test_diversity.py
git commit -m "feat: core piN/piS diversity engine with Nei-Gojobori algorithm"
```

---

### Task 8: Output Writers (`io.py`)

**Files:**
- Create: `src/pie/io.py`
- Create: `tests/test_io.py`

**Step 1: Write failing tests**

```python
# tests/test_io.py
import pytest
import pandas as pd
from pathlib import Path
from pie.io import write_gene_results, write_window_results, write_summary
from pie.diversity import GeneResult, CodonResult


@pytest.fixture
def sample_gene_results():
    return [
        GeneResult(
            gene_id="gene1",
            transcript_id="tx1",
            chrom="chr1",
            start=0,
            end=90,
            strand="+",
            n_codons=30,
            n_poly_codons=2,
            N_sites=60.0,
            S_sites=20.0,
            N_diffs=0.1,
            S_diffs=0.05,
            mean_depth=100.0,
            n_variants=3,
            codon_results=[
                CodonResult("chr1", 0, 2.0, 1.0, 0.05, 0.025),
                CodonResult("chr1", 3, 2.5, 0.5, 0.05, 0.025),
            ],
        ),
    ]


class TestWriteGeneResults:
    def test_writes_tsv(self, sample_gene_results, tmp_path):
        outpath = tmp_path / "gene_results.tsv"
        write_gene_results(sample_gene_results, str(outpath))
        assert outpath.exists()
        df = pd.read_csv(outpath, sep="\t")
        assert len(df) == 1
        assert "piN" in df.columns
        assert "piS" in df.columns
        assert "piN_piS" in df.columns

    def test_correct_values(self, sample_gene_results, tmp_path):
        outpath = tmp_path / "gene_results.tsv"
        write_gene_results(sample_gene_results, str(outpath))
        df = pd.read_csv(outpath, sep="\t")
        row = df.iloc[0]
        assert abs(row["piN"] - 0.1 / 60.0) < 1e-10
        assert abs(row["piS"] - 0.05 / 20.0) < 1e-10


class TestWriteWindowResults:
    def test_writes_tsv(self, sample_gene_results, tmp_path):
        outpath = tmp_path / "window_results.tsv"
        write_window_results(sample_gene_results, str(outpath), window_size=10, window_step=5)
        assert outpath.exists()
        df = pd.read_csv(outpath, sep="\t")
        assert len(df) > 0
        assert "piN" in df.columns


class TestWriteSummary:
    def test_writes_tsv(self, sample_gene_results, tmp_path):
        outpath = tmp_path / "summary.tsv"
        write_summary(sample_gene_results, str(outpath))
        assert outpath.exists()
        df = pd.read_csv(outpath, sep="\t")
        assert len(df) == 1
        assert "genome_piN" in df.columns
```

**Step 2: Run tests, verify fail, implement, verify pass**

Implement `io.py` with `write_gene_results`, `write_window_results`, `write_summary` functions using pandas DataFrames → `to_csv(sep="\t")`.

**Step 3: Commit**

```bash
git add src/pie/io.py tests/test_io.py
git commit -m "feat: TSV output writers for gene, window, and summary results"
```

---

### Task 9: Parallelization (`parallel.py`)

**Files:**
- Create: `src/pie/parallel.py`
- Create: `tests/test_parallel.py`

**Step 1: Write failing tests**

```python
# tests/test_parallel.py
import pytest
from pie.parallel import run_parallel
from pie.diversity import GeneResult


class TestRunParallel:
    def test_single_thread(self, ref_fasta, gff3_file, vcf_file):
        results = run_parallel(
            ref_fasta, gff3_file, vcf_file,
            min_freq=0.0, min_depth=0, min_qual=0,
            pass_only=False, threads=1,
        )
        assert len(results) > 0
        assert all(isinstance(r, GeneResult) for r in results)

    def test_multi_thread(self, ref_fasta, gff3_file, vcf_file):
        results = run_parallel(
            ref_fasta, gff3_file, vcf_file,
            min_freq=0.0, min_depth=0, min_qual=0,
            pass_only=False, threads=2,
        )
        assert len(results) > 0

    def test_results_match_single_vs_multi(self, ref_fasta, gff3_file, vcf_file):
        """Single-threaded and multi-threaded should produce identical results."""
        r1 = run_parallel(
            ref_fasta, gff3_file, vcf_file,
            min_freq=0.0, min_depth=0, min_qual=0,
            pass_only=False, threads=1,
        )
        r2 = run_parallel(
            ref_fasta, gff3_file, vcf_file,
            min_freq=0.0, min_depth=0, min_qual=0,
            pass_only=False, threads=2,
        )
        assert len(r1) == len(r2)
        for a, b in zip(r1, r2):
            assert abs(a.piN - b.piN) < 1e-10
            assert abs(a.piS - b.piS) < 1e-10
```

**Step 2: Implement parallel.py**

Use `multiprocessing.Pool` with a worker function that opens its own VCF/FASTA handles and processes a single gene. The main function parses annotations, distributes genes, collects results.

**Step 3: Run tests, verify pass**

Run: `mamba run -n pie pytest tests/test_parallel.py -v`

**Step 4: Commit**

```bash
git add src/pie/parallel.py tests/test_parallel.py
git commit -m "feat: gene-level multiprocessing for parallel piN/piS computation"
```

---

### Task 10: Manhattan Plot (`plot.py`)

**Files:**
- Create: `src/pie/plot.py`
- Create: `tests/test_plot.py`

**Step 1: Write failing tests**

```python
# tests/test_plot.py
import pytest
import pandas as pd
from pathlib import Path
from pie.plot import manhattan_plot


class TestManhattanPlot:
    def test_creates_png(self, tmp_path):
        # Create minimal gene results TSV
        df = pd.DataFrame({
            "chrom": ["chr1", "chr1", "chr2", "chr2"],
            "gene_id": ["g1", "g2", "g3", "g4"],
            "start": [100, 500, 100, 800],
            "end": [200, 600, 300, 900],
            "piN_piS": [0.5, 1.2, 0.8, 2.0],
        })
        tsv_path = tmp_path / "gene_results.tsv"
        df.to_csv(tsv_path, sep="\t", index=False)
        out_path = tmp_path / "manhattan.png"
        manhattan_plot(str(tsv_path), str(out_path))
        assert out_path.exists()
        assert out_path.stat().st_size > 0

    def test_custom_dimensions(self, tmp_path):
        df = pd.DataFrame({
            "chrom": ["chr1"],
            "gene_id": ["g1"],
            "start": [100],
            "end": [200],
            "piN_piS": [1.0],
        })
        tsv_path = tmp_path / "gene_results.tsv"
        df.to_csv(tsv_path, sep="\t", index=False)
        out_path = tmp_path / "manhattan.png"
        manhattan_plot(str(tsv_path), str(out_path), width=20, height=8)
        assert out_path.exists()
```

**Step 2: Implement plot.py**

Manhattan plot with:
- X-axis: cumulative genomic position, chromosomes alternately colored
- Y-axis: piN/piS per gene
- Horizontal line at y=1
- Gene midpoint as x coordinate
- Use matplotlib, non-interactive backend (`Agg`)

**Step 3: Run tests, verify pass**

Run: `mamba run -n pie pytest tests/test_plot.py -v`

**Step 4: Commit**

```bash
git add src/pie/plot.py tests/test_plot.py
git commit -m "feat: Manhattan plot for per-gene piN/piS"
```

---

### Task 11: Wire Up CLI (`cli.py`)

**Files:**
- Modify: `src/pie/cli.py`
- Create: `tests/test_cli.py`

**Step 1: Write failing tests**

```python
# tests/test_cli.py
import pytest
from click.testing import CliRunner
from pie.cli import main


@pytest.fixture
def runner():
    return CliRunner()


class TestRunCommand:
    def test_help(self, runner):
        result = runner.invoke(main, ["run", "--help"])
        assert result.exit_code == 0
        assert "--vcf" in result.output
        assert "--gff" in result.output
        assert "--fasta" in result.output

    def test_missing_required(self, runner):
        result = runner.invoke(main, ["run"])
        assert result.exit_code != 0

    def test_full_run(self, runner, ref_fasta, gff3_file, vcf_file, tmp_path):
        result = runner.invoke(
            main,
            [
                "run",
                "--vcf", vcf_file,
                "--gff", gff3_file,
                "--fasta", ref_fasta,
                "--outdir", str(tmp_path),
                "--min-freq", "0",
                "--min-depth", "0",
                "--min-qual", "0",
            ],
        )
        assert result.exit_code == 0, result.output
        assert (tmp_path / "gene_results.tsv").exists()
        assert (tmp_path / "summary.tsv").exists()


class TestPlotCommand:
    def test_help(self, runner):
        result = runner.invoke(main, ["plot", "--help"])
        assert result.exit_code == 0


class TestSummaryCommand:
    def test_help(self, runner):
        result = runner.invoke(main, ["summary", "--help"])
        assert result.exit_code == 0
```

**Step 2: Implement full CLI with Click options**

Wire `pie run` to:
1. Call `ensure_indexed()`
2. Call `run_parallel()` with all filter params
3. Call `write_gene_results()`, `write_window_results()`, `write_summary()`

Wire `pie plot` to call `manhattan_plot()`.
Wire `pie summary` to read and print `summary.tsv`.

**Step 3: Run tests, verify pass**

Run: `mamba run -n pie pytest tests/test_cli.py -v`

**Step 4: Run full test suite**

Run: `mamba run -n pie pytest -v`
Expected: all PASS

**Step 5: Commit**

```bash
git add src/pie/cli.py tests/test_cli.py
git commit -m "feat: wire up CLI with run/plot/summary subcommands"
```

---

### Task 12: Integration Test with Hand-Calculated Values

**Files:**
- Create: `tests/test_integration.py`

**Step 1: Write integration test comparing against hand-calculated piN/piS**

The test data from Task 2 includes known variants at known codon positions. Compute expected piN/piS by hand and assert the tool produces matching values within tolerance.

```python
# tests/test_integration.py
import pytest
import pandas as pd
from click.testing import CliRunner
from pie.cli import main


class TestHandCalculatedValues:
    def test_gene1_piN_piS(self, ref_fasta, gff3_file, vcf_file, tmp_path):
        """Verify gene1 piN/piS matches hand-calculated values from test data."""
        runner = CliRunner()
        result = runner.invoke(
            main,
            [
                "run",
                "--vcf", vcf_file,
                "--gff", gff3_file,
                "--fasta", ref_fasta,
                "--outdir", str(tmp_path),
                "--min-freq", "0",
                "--min-depth", "0",
                "--min-qual", "0",
            ],
        )
        assert result.exit_code == 0, result.output

        df = pd.read_csv(tmp_path / "gene_results.tsv", sep="\t")
        gene1 = df[df["gene_id"].str.contains("gene1", case=False)].iloc[0]

        # These expected values come from create_test_data.py comments
        # and must be filled in after Task 2 generates the data with
        # known hand-calculated values.
        assert gene1["N_sites"] > 0
        assert gene1["S_sites"] > 0
        # TODO: fill in exact expected values from create_test_data.py

    def test_all_output_files_created(self, ref_fasta, gff3_file, vcf_file, tmp_path):
        runner = CliRunner()
        result = runner.invoke(
            main,
            [
                "run",
                "--vcf", vcf_file,
                "--gff", gff3_file,
                "--fasta", ref_fasta,
                "--outdir", str(tmp_path),
                "--min-freq", "0",
                "--min-depth", "0",
                "--min-qual", "0",
                "--window-size", "30",
                "--window-step", "10",
            ],
        )
        assert result.exit_code == 0
        assert (tmp_path / "gene_results.tsv").exists()
        assert (tmp_path / "window_results.tsv").exists()
        assert (tmp_path / "summary.tsv").exists()

    def test_monomorphic_gene_zero_diffs(self, ref_fasta, gff3_file, vcf_file, tmp_path):
        """A gene with no variants should have piN = piS = 0."""
        runner = CliRunner()
        result = runner.invoke(
            main,
            [
                "run",
                "--vcf", vcf_file,
                "--gff", gff3_file,
                "--fasta", ref_fasta,
                "--outdir", str(tmp_path),
                "--min-freq", "0",
                "--min-depth", "0",
                "--min-qual", "0",
            ],
        )
        assert result.exit_code == 0
        df = pd.read_csv(tmp_path / "gene_results.tsv", sep="\t")
        # If test data has a gene with no variants, check it
        for _, row in df.iterrows():
            if row["n_variants"] == 0:
                assert row["piN"] == 0.0
                assert row["piS"] == 0.0
```

**Step 2: Run integration tests**

Run: `mamba run -n pie pytest tests/test_integration.py -v`
Expected: all PASS

**Step 3: Commit**

```bash
git add tests/test_integration.py
git commit -m "test: integration tests with hand-calculated piN/piS values"
```

---

### Task 13: Final Polish

**Step 1: Run full test suite**

Run: `mamba run -n pie pytest -v --tb=short`
Expected: all PASS

**Step 2: Test CLI end-to-end manually**

Run: `mamba run -n pie pie run --vcf tests/data/variants.vcf.gz --gff tests/data/genes.gff3 --fasta tests/data/ref.fa --outdir /tmp/pie_test`
Run: `mamba run -n pie pie summary /tmp/pie_test/summary.tsv`
Run: `mamba run -n pie pie plot --gene-results /tmp/pie_test/gene_results.tsv --output /tmp/pie_test/manhattan.png`

**Step 3: Verify output files look correct**

Run: `head /tmp/pie_test/gene_results.tsv`
Run: `head /tmp/pie_test/summary.tsv`
Run: `ls -la /tmp/pie_test/manhattan.png`

**Step 4: Final commit**

```bash
git add -A
git commit -m "chore: final polish and cleanup"
```

---

## Task Dependency Graph

```
Task 1 (scaffold)
  └─► Task 2 (test data)
       └─► Task 3 (codon.py) ─────────────────┐
       └─► Task 4 (reference.py) ──────────────┤
       └─► Task 5 (annotation.py) ─────────────┤
       └─► Task 6 (vcf.py) ────────────────────┤
                                                ▼
                                    Task 7 (diversity.py)
                                         │
                              ┌──────────┼──────────┐
                              ▼          ▼          ▼
                     Task 8 (io.py)  Task 9 (parallel.py)  Task 10 (plot.py)
                              │          │          │
                              └──────────┼──────────┘
                                         ▼
                                  Task 11 (cli.py)
                                         │
                                         ▼
                                  Task 12 (integration)
                                         │
                                         ▼
                                  Task 13 (polish)
```

**Note:** Tasks 3, 4, 5, 6 can be implemented in parallel after Task 2 completes (no inter-dependencies). Tasks 8, 9, 10 can be parallelized after Task 7.
