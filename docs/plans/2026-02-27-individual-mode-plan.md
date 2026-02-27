# Individual Sequencing Mode Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add `--mode individual` support to `pie run`, deriving pooled allele frequencies from GT fields across multiple diploid samples.

**Architecture:** New `IndividualVariantReader` class in `vcf.py` with the same `fetch() -> list[Variant]` interface as `VariantReader`. Downstream code (`diversity.py`, `parallel.py`) is unchanged. CLI adds `--mode` flag with validation. Output writers conditionally include sample metadata columns.

**Tech Stack:** Python 3.12+, Click CLI, cyvcf2, pytest

**Design doc:** `docs/plans/2026-02-27-individual-mode-design.md`

---

### Task 1: Add `call_rate` field to `Variant` dataclass

Extend `Variant` to carry per-site call rate (set in individual mode, `None` in pool mode).

**Files:**
- Modify: `src/pie/vcf.py:13-19`
- Test: `tests/test_vcf.py`

**Step 1: Add optional `call_rate` field**

In `src/pie/vcf.py`, add `call_rate` to the `Variant` dataclass:

```python
@dataclass(slots=True)
class Variant:
    pos: int      # 0-based genomic position
    ref: str      # reference allele
    alt: str      # alternate allele
    freq: float   # alt allele frequency (0-1)
    depth: int    # total depth (read depth in pool mode, AN in individual mode)
    call_rate: float | None = None  # fraction of called samples (individual mode only)
```

**Step 2: Run existing tests to verify backward compatibility**

Run: `mamba run -n pie pytest tests/test_vcf.py tests/test_integration.py -v`
Expected: All PASS (default `None` is backward-compatible).

**Step 3: Commit**

```bash
git add src/pie/vcf.py
git commit -m "feat: add optional call_rate field to Variant dataclass"
```

---

### Task 2: Add `n_samples` and `call_rates` to `GeneResult`; collect in `compute_gene_diversity`

**Files:**
- Modify: `src/pie/diversity.py:50-68` (GeneResult), `src/pie/diversity.py:241-357` (compute_gene_diversity)
- Test: `tests/test_diversity.py`, `tests/test_integration.py`

**Step 1: Add optional fields to `GeneResult`**

In `src/pie/diversity.py`, add fields after `codon_results`:

```python
@dataclass
class GeneResult:
    """Per-gene diversity result."""

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
    mean_variant_depth: float
    n_variants: int
    codon_results: list[CodonResult] = field(default_factory=list)
    n_samples: int | None = None
    call_rates: list[float] | None = None

    @property
    def piN(self) -> float:
        return self.N_diffs / self.N_sites if self.N_sites > 0 else 0.0

    @property
    def piS(self) -> float:
        return self.S_diffs / self.S_sites if self.S_sites > 0 else 0.0

    @property
    def piN_piS(self) -> float | None:
        return self.piN / self.piS if self.piS > 0 else None

    @property
    def mean_call_rate(self) -> float | None:
        if self.call_rates is None or len(self.call_rates) == 0:
            return None
        return sum(self.call_rates) / len(self.call_rates)
```

**Step 2: Collect call_rates in `compute_gene_diversity`**

At the end of `compute_gene_diversity`, before building the return value, collect call rates from variants:

```python
    # Collect per-variant call rates (individual mode only)
    cr_list = [v.call_rate for v in all_variants if v.call_rate is not None]

    return GeneResult(
        ...
        codon_results=codon_results,
        call_rates=cr_list if cr_list else None,
    )
```

Note: `n_samples` is NOT set here — it's set by the parallel runner (which knows the reader's sample count). This keeps `compute_gene_diversity` agnostic to the reader type.

**Step 3: Run existing tests**

Run: `mamba run -n pie pytest tests/test_diversity.py tests/test_integration.py -v`
Expected: All PASS.

**Step 4: Commit**

```bash
git add src/pie/diversity.py
git commit -m "feat: add n_samples, call_rates, mean_call_rate to GeneResult"
```

---

### Task 3: Create multi-sample individual VCF test fixture

**Files:**
- Modify: `tests/conftest.py`

**Step 1: Add `individual_vcf_file` fixture**

4 diploid samples (S1–S4), same variant positions as existing test data, with varying missingness:

```python
@pytest.fixture
def individual_vcf_file(tmp_path):
    """Multi-sample VCF with GT fields for individual-mode testing.

    4 diploid samples (S1-S4), same positions as test data:
      pos 6   T>C  S1:0/1  S2:0/0  S3:0/1  S4:./.  → called=3, AN=6, AC=2, freq=1/3, call_rate=0.75
      pos 7   G>A  S1:0/0  S2:0/1  S3:0/1  S4:0/0  → called=4, AN=8, AC=2, freq=1/4, call_rate=1.00
      pos 195 A>T  S1:0/1  S2:0/1  S3:1/1  S4:0/1  → called=4, AN=8, AC=6, freq=3/4, call_rate=1.00
      pos 297 A>G  S1:0/1  S2:./.  S3:./.  S4:0/0  → called=2, AN=4, AC=1, freq=1/4, call_rate=0.50
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
```

**Step 2: Add `individual_multiallelic_vcf_file` fixture**

Multi-sample VCF with a non-decomposed multiallelic site:

```python
@pytest.fixture
def individual_multiallelic_vcf_file(tmp_path):
    """Multi-sample VCF with a multiallelic site for individual-mode testing.

      pos 6   T>C    S1:0/1  S2:0/0  S3:0/1  S4:0/0  → AN=8, AC=2, freq=1/4
      pos 7   G>A,C  S1:0/1  S2:0/2  S3:0/1  S4:1/2  → AN=8, AC_A=3, AC_C=3, freq_A=3/8, freq_C=3/8
      pos 195 A>T    S1:0/1  S2:0/0  S3:0/1  S4:0/1  → AN=8, AC=3, freq=3/8
    """
    vcf_content = """\
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=chr1,length=350>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4
chr1\t6\t.\tT\tC\t30\t.\t.\tGT\t0/1\t0/0\t0/1\t0/0
chr1\t7\t.\tG\tA,C\t50\t.\t.\tGT\t0/1\t0/2\t0/1\t1/2
chr1\t195\t.\tA\tT\t45\t.\t.\tGT\t0/1\t0/0\t0/1\t0/1
"""
    vcf_path = tmp_path / "individual_multiallelic.vcf"
    vcf_path.write_text(vcf_content)
    return _bgzip_and_index(vcf_path)
```

**Step 3: Commit**

```bash
git add tests/conftest.py
git commit -m "test: add multi-sample individual VCF fixtures"
```

---

### Task 4: `IndividualVariantReader` — basic GT→frequency extraction

**Files:**
- Modify: `src/pie/vcf.py`
- Test: `tests/test_vcf.py`

**Step 1: Write failing tests**

Add a new test class `TestIndividualVariantReader` in `tests/test_vcf.py`:

```python
from pie.vcf import IndividualVariantReader


class TestIndividualVariantReader:
    def test_basic_gt_frequency(self, individual_vcf_file):
        """All 4 samples, no filtering. Verify GT-derived frequencies."""
        with IndividualVariantReader(
            individual_vcf_file,
            samples=["S1", "S2", "S3", "S4"],
            min_freq=0.0, min_qual=0.0, pass_only=False,
            keep_multiallelic=False, min_call_rate=0.0, min_an=0,
        ) as reader:
            variants = reader.fetch("chr1", 0, 350)
            # 4 biallelic sites, all pass with no filters
            assert len(variants) == 4

            v6 = next(v for v in variants if v.pos == 5)
            assert v6.ref == "T" and v6.alt == "C"
            assert abs(v6.freq - 2 / 6) < 1e-6   # AC=2, AN=6
            assert v6.depth == 6                    # AN
            assert abs(v6.call_rate - 0.75) < 1e-6  # 3/4

            v7 = next(v for v in variants if v.pos == 6)
            assert abs(v7.freq - 2 / 8) < 1e-6
            assert v7.depth == 8
            assert abs(v7.call_rate - 1.0) < 1e-6

            v195 = next(v for v in variants if v.pos == 194)
            assert abs(v195.freq - 6 / 8) < 1e-6  # 3 hets + 1 hom_alt = 6 alt alleles
            assert v195.depth == 8

            v297 = next(v for v in variants if v.pos == 296)
            assert abs(v297.freq - 1 / 4) < 1e-6  # 1 het in 2 called = AC=1, AN=4
            assert v297.depth == 4
            assert abs(v297.call_rate - 0.50) < 1e-6

    def test_context_manager(self, individual_vcf_file):
        with IndividualVariantReader(
            individual_vcf_file, samples=["S1", "S2", "S3", "S4"],
            min_freq=0.0, min_qual=0.0, pass_only=False,
            keep_multiallelic=False, min_call_rate=0.0, min_an=0,
        ) as reader:
            assert reader is not None
            assert reader.n_samples == 4
```

**Step 2: Run tests to verify they fail**

Run: `mamba run -n pie pytest tests/test_vcf.py::TestIndividualVariantReader -v`
Expected: FAIL (ImportError — `IndividualVariantReader` does not exist).

**Step 3: Implement `IndividualVariantReader`**

Add to `src/pie/vcf.py`:

```python
class IndividualVariantReader:
    """Variant reader for individual-sequencing data (GT-based frequencies).

    Derives pooled allele frequencies from genotype fields across selected
    diploid samples.  Exposes the same ``fetch()`` interface as
    ``VariantReader`` so downstream code is unchanged.
    """

    def __init__(self, vcf_path: str, samples: list[str],
                 min_freq: float = 0.01, min_qual: float = 20.0,
                 pass_only: bool = False, keep_multiallelic: bool = False,
                 min_call_rate: float = 0.8, min_an: int = 2):
        self._min_freq = min_freq
        self._min_qual = min_qual
        self._pass_only = pass_only
        self._keep_multiallelic = keep_multiallelic
        self._min_call_rate = min_call_rate
        self._min_an = min_an
        self._vcf = VCF(vcf_path, samples=samples)
        self._n_samples = len(samples)

    def close(self):
        self._vcf.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    @property
    def n_samples(self) -> int:
        return self._n_samples

    def fetch(self, chrom: str, start: int, end: int) -> list[Variant]:
        """Fetch filtered variants in region (0-based, half-open).

        Per-record logic (single GT pass):
        1. QUAL/PASS/REF filters
        2. One pass through selected samples' GT:
           - Count called samples, accumulate ref/alt allele counts
           - AN = 2 * called (diploid)
        3. Check call_rate >= min_call_rate and AN >= min_an
        4. For each ALT with alt_count > 0: freq = alt_count / AN
        5. Multiallelic grouping, then min_freq filter
        """
        region = f"{chrom}:{start + 1}-{end}"

        # (pos, ref, alt, freq, AN, ref_count, alt_count, call_rate)
        raw: list[tuple[int, str, str, float, int, int, int, float]] = []
        multiallelic_pos: set[int] = set()
        seen_pos: dict[int, int] = {}

        for record in self._vcf(region):
            if self._pass_only and record.FILTER is not None:
                continue
            if record.QUAL is not None and record.QUAL < self._min_qual:
                continue
            if record.REF not in _VALID_BASES:
                continue

            pos0 = record.POS - 1

            # Track multiallelic positions
            n_valid = sum(1 for a in record.ALT if a in _VALID_BASES)
            seen_pos[pos0] = seen_pos.get(pos0, 0) + n_valid
            if seen_pos[pos0] > 1:
                multiallelic_pos.add(pos0)

            # --- Single-pass GT extraction ---
            n_alts = len(record.ALT)
            called = 0
            ref_count = 0
            alt_counts = [0] * n_alts

            for gt in record.genotypes:  # [allele1, allele2, is_phased]
                a1, a2 = gt[0], gt[1]
                if a1 < 0 or a2 < 0:
                    continue
                called += 1
                for allele in (a1, a2):
                    if allele == 0:
                        ref_count += 1
                    elif 1 <= allele <= n_alts:
                        alt_counts[allele - 1] += 1

            if called == 0:
                continue

            call_rate = called / self._n_samples
            if call_rate < self._min_call_rate:
                continue

            an = 2 * called
            if an < self._min_an:
                continue

            for alt_idx, alt_allele in enumerate(record.ALT):
                if alt_allele not in _VALID_BASES:
                    continue
                ac = alt_counts[alt_idx]
                if ac == 0:
                    continue
                freq = ac / an
                raw.append((pos0, record.REF, alt_allele, freq, an,
                            ref_count, ac, call_rate))

        if not raw:
            return []

        # Group by position for multiallelic handling
        by_pos: dict[int, list[tuple]] = {}
        for rec in raw:
            by_pos.setdefault(rec[0], []).append(rec)

        variants: list[Variant] = []
        for pos in sorted(by_pos):
            group = by_pos[pos]
            is_multiallelic = pos in multiallelic_pos
            if not is_multiallelic:
                p, ref, alt, freq, depth, _, _, cr = group[0]
                if freq < self._min_freq:
                    continue
                variants.append(Variant(pos=p, ref=ref, alt=alt,
                                        freq=freq, depth=depth, call_rate=cr))
            elif not self._keep_multiallelic:
                continue
            else:
                for r in group:
                    p, ref, alt, freq, depth, _, _, cr = r
                    if freq >= self._min_freq:
                        variants.append(Variant(pos=p, ref=ref, alt=alt,
                                                freq=freq, depth=depth,
                                                call_rate=cr))

        return variants
```

**Step 4: Run tests to verify they pass**

Run: `mamba run -n pie pytest tests/test_vcf.py::TestIndividualVariantReader -v`
Expected: All PASS.

**Step 5: Commit**

```bash
git add src/pie/vcf.py tests/test_vcf.py
git commit -m "feat: add IndividualVariantReader with GT-based frequency extraction"
```

---

### Task 5: `IndividualVariantReader` — filtering tests

**Files:**
- Test: `tests/test_vcf.py`

**Step 1: Write and run filtering tests**

Add to `TestIndividualVariantReader` in `tests/test_vcf.py`:

```python
    def test_min_call_rate_filter(self, individual_vcf_file):
        """min_call_rate=0.8 filters pos 6 (0.75) and pos 297 (0.50)."""
        with IndividualVariantReader(
            individual_vcf_file, samples=["S1", "S2", "S3", "S4"],
            min_freq=0.0, min_qual=0.0, pass_only=False,
            keep_multiallelic=False, min_call_rate=0.8, min_an=0,
        ) as reader:
            variants = reader.fetch("chr1", 0, 350)
            positions = [v.pos for v in variants]
            assert 5 not in positions   # pos 6: call_rate=0.75 < 0.8
            assert 6 in positions       # pos 7: call_rate=1.0
            assert 194 in positions     # pos 195: call_rate=1.0
            assert 296 not in positions  # pos 297: call_rate=0.50 < 0.8
            assert len(variants) == 2

    def test_min_an_filter(self, individual_vcf_file):
        """min_an=6 filters pos 297 (AN=4) but keeps pos 6 (AN=6)."""
        with IndividualVariantReader(
            individual_vcf_file, samples=["S1", "S2", "S3", "S4"],
            min_freq=0.0, min_qual=0.0, pass_only=False,
            keep_multiallelic=False, min_call_rate=0.0, min_an=6,
        ) as reader:
            variants = reader.fetch("chr1", 0, 350)
            positions = [v.pos for v in variants]
            assert 5 in positions       # AN=6 >= 6
            assert 6 in positions       # AN=8 >= 6
            assert 194 in positions     # AN=8 >= 6
            assert 296 not in positions  # AN=4 < 6

    def test_min_freq_on_gt_derived_af(self, individual_vcf_file):
        """min_freq=0.30 filters pos 7 (freq=0.25) and pos 297 (freq=0.25)."""
        with IndividualVariantReader(
            individual_vcf_file, samples=["S1", "S2", "S3", "S4"],
            min_freq=0.30, min_qual=0.0, pass_only=False,
            keep_multiallelic=False, min_call_rate=0.0, min_an=0,
        ) as reader:
            variants = reader.fetch("chr1", 0, 350)
            positions = [v.pos for v in variants]
            assert 5 in positions       # freq=1/3=0.333 >= 0.30
            assert 6 not in positions   # freq=1/4=0.25 < 0.30
            assert 194 in positions     # freq=3/4=0.75 >= 0.30
            assert 296 not in positions  # freq=1/4=0.25 < 0.30

    def test_qual_filter(self, individual_vcf_file):
        """min_qual=20 filters pos 297 (QUAL=15)."""
        with IndividualVariantReader(
            individual_vcf_file, samples=["S1", "S2", "S3", "S4"],
            min_freq=0.0, min_qual=20.0, pass_only=False,
            keep_multiallelic=False, min_call_rate=0.0, min_an=0,
        ) as reader:
            variants = reader.fetch("chr1", 0, 350)
            positions = [v.pos for v in variants]
            assert 296 not in positions  # QUAL=15 < 20

    def test_sample_subset(self, individual_vcf_file):
        """Selecting only S1,S2: pos 6 T>C has AC=1/AN=4, call_rate=1.0."""
        with IndividualVariantReader(
            individual_vcf_file, samples=["S1", "S2"],
            min_freq=0.0, min_qual=0.0, pass_only=False,
            keep_multiallelic=False, min_call_rate=0.0, min_an=0,
        ) as reader:
            assert reader.n_samples == 2
            variants = reader.fetch("chr1", 0, 100)
            v6 = next(v for v in variants if v.pos == 5)
            # S1: 0/1 (AC=1), S2: 0/0 (AC=0) → AC=1, AN=4
            assert abs(v6.freq - 1 / 4) < 1e-6
            assert v6.depth == 4
            assert abs(v6.call_rate - 1.0) < 1e-6  # both called

    def test_all_missing_site_skipped(self, individual_vcf_file):
        """If all selected samples are missing, the site is skipped."""
        # S4 is ./. at pos 6 — select only S4
        with IndividualVariantReader(
            individual_vcf_file, samples=["S4"],
            min_freq=0.0, min_qual=0.0, pass_only=False,
            keep_multiallelic=False, min_call_rate=0.0, min_an=0,
        ) as reader:
            variants = reader.fetch("chr1", 0, 10)
            # pos 6: S4 is ./., so 0 called → skipped
            positions = [v.pos for v in variants]
            assert 5 not in positions
```

**Step 2: Run tests**

Run: `mamba run -n pie pytest tests/test_vcf.py::TestIndividualVariantReader -v`
Expected: All PASS (implementation from Task 4 already handles these).

**Step 3: Commit**

```bash
git add tests/test_vcf.py
git commit -m "test: add filtering and sample subset tests for IndividualVariantReader"
```

---

### Task 6: `IndividualVariantReader` — multiallelic handling

**Files:**
- Test: `tests/test_vcf.py`

**Step 1: Write and run multiallelic tests**

Add to `tests/test_vcf.py`:

```python
class TestIndividualMultiallelic:
    def test_default_skips_multiallelic(self, individual_multiallelic_vcf_file):
        """pos 7 has ALT=A,C → multiallelic, skipped by default."""
        with IndividualVariantReader(
            individual_multiallelic_vcf_file,
            samples=["S1", "S2", "S3", "S4"],
            min_freq=0.0, min_qual=0.0, pass_only=False,
            keep_multiallelic=False, min_call_rate=0.0, min_an=0,
        ) as reader:
            variants = reader.fetch("chr1", 0, 350)
            positions = [v.pos for v in variants]
            assert 6 not in positions  # multiallelic skipped
            assert 5 in positions
            assert 194 in positions
            assert len(variants) == 2

    def test_keep_multiallelic(self, individual_multiallelic_vcf_file):
        """With keep_multiallelic=True, both ALTs at pos 7 are kept.

        pos 7 G>A,C: S1:0/1 S2:0/2 S3:0/1 S4:1/2
        AN=8, AC_A=3, AC_C=3, freq_A=3/8=0.375, freq_C=3/8=0.375
        """
        with IndividualVariantReader(
            individual_multiallelic_vcf_file,
            samples=["S1", "S2", "S3", "S4"],
            min_freq=0.0, min_qual=0.0, pass_only=False,
            keep_multiallelic=True, min_call_rate=0.0, min_an=0,
        ) as reader:
            variants = reader.fetch("chr1", 0, 350)
            ma = [v for v in variants if v.pos == 6]
            assert len(ma) == 2
            va = next(v for v in ma if v.alt == "A")
            vc = next(v for v in ma if v.alt == "C")
            assert abs(va.freq - 3 / 8) < 1e-6
            assert abs(vc.freq - 3 / 8) < 1e-6
            assert va.depth == 8  # AN
            assert vc.depth == 8

    def test_multiallelic_min_freq_filter(self, individual_multiallelic_vcf_file):
        """min_freq filters individual ALTs at multiallelic site."""
        with IndividualVariantReader(
            individual_multiallelic_vcf_file,
            samples=["S1", "S2", "S3", "S4"],
            min_freq=0.40, min_qual=0.0, pass_only=False,
            keep_multiallelic=True, min_call_rate=0.0, min_an=0,
        ) as reader:
            variants = reader.fetch("chr1", 0, 350)
            # Both ALTs at pos 7 have freq=0.375 < 0.40
            ma = [v for v in variants if v.pos == 6]
            assert len(ma) == 0
```

**Step 2: Run tests**

Run: `mamba run -n pie pytest tests/test_vcf.py::TestIndividualMultiallelic -v`
Expected: All PASS.

**Step 3: Commit**

```bash
git add tests/test_vcf.py
git commit -m "test: add multiallelic tests for IndividualVariantReader"
```

---

### Task 7: Modify output writers for conditional metadata columns

**Files:**
- Modify: `src/pie/io.py`
- Test: `tests/test_io.py`

**Step 1: Write failing tests**

Add to `tests/test_io.py`:

```python
class TestIndividualModeColumns:
    def _make_result(self, n_samples=None, call_rates=None):
        """Helper to create a GeneResult with optional ind-mode metadata."""
        return GeneResult(
            gene_id="gene1", transcript_id="tx1", chrom="chr1",
            start=0, end=90, strand="+", n_codons=29, n_poly_codons=2,
            N_sites=59.0, S_sites=27.0, N_diffs=0.42, S_diffs=0.32,
            mean_variant_depth=100.0, n_variants=2,
            n_samples=n_samples, call_rates=call_rates,
        )

    def test_pool_mode_no_extra_columns(self, tmp_path):
        """Pool mode: n_samples=None → no n_samples/mean_call_rate columns."""
        results = [self._make_result()]
        path = str(tmp_path / "gene_results.tsv")
        write_gene_results(results, path)
        df = pd.read_csv(path, sep="\t")
        assert "n_samples" not in df.columns
        assert "mean_call_rate" not in df.columns

    def test_individual_mode_extra_columns(self, tmp_path):
        """Individual mode: n_samples set → extra columns present."""
        results = [self._make_result(n_samples=4, call_rates=[0.75, 1.0])]
        path = str(tmp_path / "gene_results.tsv")
        write_gene_results(results, path)
        df = pd.read_csv(path, sep="\t")
        assert "n_samples" in df.columns
        assert "mean_call_rate" in df.columns
        assert df.iloc[0]["n_samples"] == 4
        assert abs(df.iloc[0]["mean_call_rate"] - 0.875) < 1e-6  # (0.75+1.0)/2

    def test_summary_pool_mode(self, tmp_path):
        """Pool mode summary: no n_samples_selected/mean_call_rate."""
        results = [self._make_result()]
        path = str(tmp_path / "summary.tsv")
        write_summary(results, path)
        df = pd.read_csv(path, sep="\t")
        assert "n_samples_selected" not in df.columns
        assert "mean_call_rate" not in df.columns

    def test_summary_individual_mode(self, tmp_path):
        """Individual mode summary: includes n_samples_selected and mean_call_rate."""
        r1 = self._make_result(n_samples=4, call_rates=[0.75, 1.0])
        r2 = self._make_result(n_samples=4, call_rates=[1.0])
        results = [r1, r2]
        path = str(tmp_path / "summary.tsv")
        write_summary(results, path)
        df = pd.read_csv(path, sep="\t")
        assert df.iloc[0]["n_samples_selected"] == 4
        # Variant-site-weighted: (0.75 + 1.0 + 1.0) / 3 = 0.9167
        assert abs(df.iloc[0]["mean_call_rate"] - (0.75 + 1.0 + 1.0) / 3) < 1e-4
```

**Step 2: Run tests to verify they fail**

Run: `mamba run -n pie pytest tests/test_io.py::TestIndividualModeColumns -v`
Expected: FAIL (columns not yet added).

**Step 3: Implement conditional columns in `write_gene_results`**

In `src/pie/io.py`, modify `write_gene_results`:

```python
def write_gene_results(results: list[GeneResult], path: str) -> None:
    """Write per-gene TSV with piN, piS, and piN/piS columns."""
    ind_mode = any(r.n_samples is not None for r in results)
    rows = []
    for r in results:
        piN = r.piN
        piS = r.piS
        piN_piS = r.piN_piS
        row = {
            "chrom": r.chrom,
            "gene_id": r.gene_id,
            "transcript_id": r.transcript_id,
            "start": r.start,
            "end": r.end,
            "strand": r.strand,
            "n_codons": r.n_codons,
            "n_poly_codons": r.n_poly_codons,
            "N_sites": r.N_sites,
            "S_sites": r.S_sites,
            "N_diffs": r.N_diffs,
            "S_diffs": r.S_diffs,
            "piN": piN,
            "piS": piS,
            "piN_piS": piN_piS if piN_piS is not None else "NA",
            "mean_variant_depth": r.mean_variant_depth,
            "n_variants": r.n_variants,
        }
        if ind_mode:
            row["n_samples"] = r.n_samples
            row["mean_call_rate"] = r.mean_call_rate
        rows.append(row)
    df = pd.DataFrame(rows)
    if not df.empty:
        df["piN_piS"] = pd.to_numeric(df["piN_piS"], errors="coerce")
    df.to_csv(path, sep="\t", index=False)
```

**Step 4: Implement conditional columns in `write_summary`**

In `src/pie/io.py`, modify `write_summary` — add after the existing `row` dict:

```python
def write_summary(results: list[GeneResult], path: str) -> None:
    """Write single-row genome-wide summary TSV."""
    total_genes = len(results)
    total_codons = sum(r.n_codons for r in results)
    total_variants = sum(r.n_variants for r in results)

    total_N_sites = sum(r.N_sites for r in results)
    total_S_sites = sum(r.S_sites for r in results)
    total_N_diffs = sum(r.N_diffs for r in results)
    total_S_diffs = sum(r.S_diffs for r in results)

    genome_piN = total_N_diffs / total_N_sites if total_N_sites > 0 else 0.0
    genome_piS = total_S_diffs / total_S_sites if total_S_sites > 0 else 0.0
    genome_piN_piS = genome_piN / genome_piS if genome_piS > 0 else None

    gene_piNs = [r.piN for r in results]
    gene_piSs = [r.piS for r in results]

    row = {
        "total_genes": total_genes,
        "total_codons": total_codons,
        "total_variants": total_variants,
        "genome_piN": genome_piN,
        "genome_piS": genome_piS,
        "genome_piN_piS": genome_piN_piS if genome_piN_piS is not None else "NA",
        "mean_gene_piN": float(np.mean(gene_piNs)) if gene_piNs else 0.0,
        "mean_gene_piS": float(np.mean(gene_piSs)) if gene_piSs else 0.0,
        "median_gene_piN": float(np.median(gene_piNs)) if gene_piNs else 0.0,
        "median_gene_piS": float(np.median(gene_piSs)) if gene_piSs else 0.0,
    }

    # Individual mode metadata
    ind_mode = any(r.n_samples is not None for r in results)
    if ind_mode:
        row["n_samples_selected"] = results[0].n_samples if results else 0
        # Variant-site-weighted mean call rate
        all_call_rates = []
        for r in results:
            if r.call_rates:
                all_call_rates.extend(r.call_rates)
        row["mean_call_rate"] = (
            sum(all_call_rates) / len(all_call_rates) if all_call_rates else 0.0
        )

    df = pd.DataFrame([row])
    df["genome_piN_piS"] = pd.to_numeric(df["genome_piN_piS"], errors="coerce")
    df.to_csv(path, sep="\t", index=False)
```

**Step 5: Run tests**

Run: `mamba run -n pie pytest tests/test_io.py -v`
Expected: All PASS (including existing pool-mode tests).

**Step 6: Commit**

```bash
git add src/pie/io.py tests/test_io.py
git commit -m "feat: add conditional n_samples/mean_call_rate columns to output writers"
```

---

### Task 8: Wire parallel runner for mode selection

**Files:**
- Modify: `src/pie/parallel.py`
- Test: `tests/test_parallel.py`

**Step 1: Read current `tests/test_parallel.py`**

Check existing tests so we don't break them.

**Step 2: Modify `run_parallel` and `_worker_init`**

In `src/pie/parallel.py`:

```python
"""Gene-level multiprocessing for parallel piN/piS computation."""

import logging
from multiprocessing import Pool

from pie.annotation import GeneModel, parse_annotations
from pie.reference import ReferenceGenome
from pie.vcf import VariantReader, IndividualVariantReader
from pie.diversity import compute_gene_diversity, GeneResult

log = logging.getLogger(__name__)


def _worker_init(fasta_path, vcf_path, min_freq, min_depth, min_qual,
                 pass_only, keep_multiallelic, exclude_stops, sample,
                 mode, samples, min_call_rate, min_an):
    """Initialize per-worker file handles (stored in globals)."""
    global _ref, _vcf, _exclude_stops, _n_samples
    _ref = ReferenceGenome(fasta_path)
    _exclude_stops = exclude_stops

    if mode == "individual":
        _vcf = IndividualVariantReader(
            vcf_path, samples=samples, min_freq=min_freq,
            min_qual=min_qual, pass_only=pass_only,
            keep_multiallelic=keep_multiallelic,
            min_call_rate=min_call_rate, min_an=min_an,
        )
        _n_samples = _vcf.n_samples
    else:
        _vcf = VariantReader(
            vcf_path, min_freq=min_freq, min_depth=min_depth,
            min_qual=min_qual, pass_only=pass_only,
            keep_multiallelic=keep_multiallelic, sample=sample,
        )
        _n_samples = None


def _worker_cleanup():
    """Close per-worker file handles."""
    global _ref, _vcf
    _ref.close()
    _vcf.close()


def _process_gene(gene: GeneModel) -> GeneResult:
    """Process a single gene using worker-local handles."""
    result = compute_gene_diversity(gene, _ref, _vcf, exclude_stops=_exclude_stops)
    result.n_samples = _n_samples
    return result


def run_parallel(
    fasta_path: str,
    gff_path: str,
    vcf_path: str,
    min_freq: float = 0.01,
    min_depth: int = 10,
    min_qual: float = 20.0,
    pass_only: bool = False,
    keep_multiallelic: bool = False,
    exclude_stops: bool = True,
    threads: int = 1,
    sample: str | None = None,
    mode: str = "pool",
    samples: list[str] | None = None,
    min_call_rate: float = 0.8,
    min_an: int = 2,
) -> list[GeneResult]:
    """Run piN/piS analysis across all genes.

    Returns list of GeneResult sorted by (chrom, start).
    """
    genes = parse_annotations(gff_path)
    log.info("Parsed %d genes from %s", len(genes), gff_path)

    init_args = (fasta_path, vcf_path, min_freq, min_depth, min_qual,
                 pass_only, keep_multiallelic, exclude_stops, sample,
                 mode, samples, min_call_rate, min_an)

    if threads <= 1:
        _worker_init(*init_args)
        try:
            results = [_process_gene(g) for g in genes]
        finally:
            _worker_cleanup()
    else:
        with Pool(
            processes=threads,
            initializer=_worker_init,
            initargs=init_args,
        ) as pool:
            results = pool.map(_process_gene, genes)

    results.sort(key=lambda r: (r.chrom, r.start))
    return results
```

**Step 3: Run existing tests**

Run: `mamba run -n pie pytest tests/test_parallel.py tests/test_integration.py -v`
Expected: All PASS (default `mode="pool"` preserves existing behavior).

**Step 4: Commit**

```bash
git add src/pie/parallel.py
git commit -m "feat: wire parallel runner for pool/individual mode selection"
```

---

### Task 9: CLI `--mode` option and validation

**Files:**
- Modify: `src/pie/cli.py`
- Test: `tests/test_cli.py`

**Step 1: Write failing CLI validation tests**

Add to `tests/test_cli.py`:

```python
class TestModeValidation:
    def test_help_shows_mode(self, runner):
        result = runner.invoke(main, ["run", "--help"])
        assert "--mode" in result.output

    def test_default_mode_is_pool(self, runner, ref_fasta, gff3_file, vcf_file, tmp_path):
        """No --mode → defaults to pool, works as before."""
        result = runner.invoke(main, [
            "run", "--vcf", vcf_file, "--gff", gff3_file,
            "--fasta", ref_fasta, "--outdir", str(tmp_path),
            "--min-freq", "0", "--min-depth", "0", "--min-qual", "0",
        ])
        assert result.exit_code == 0, result.output

    def test_mode_pool_explicit(self, runner, ref_fasta, gff3_file, vcf_file, tmp_path):
        result = runner.invoke(main, [
            "run", "--mode", "pool", "--vcf", vcf_file, "--gff", gff3_file,
            "--fasta", ref_fasta, "--outdir", str(tmp_path),
            "--min-freq", "0", "--min-depth", "0", "--min-qual", "0",
        ])
        assert result.exit_code == 0, result.output

    def test_mode_ind_alias(self, runner, ref_fasta, gff3_file, individual_vcf_file, tmp_path):
        """'ind' is accepted as alias for 'individual'."""
        result = runner.invoke(main, [
            "run", "--mode", "ind", "--vcf", individual_vcf_file,
            "--gff", gff3_file, "--fasta", ref_fasta,
            "--outdir", str(tmp_path),
            "--min-freq", "0", "--min-qual", "0", "--min-call-rate", "0",
        ])
        assert result.exit_code == 0, result.output

    def test_pool_with_samples_error(self, runner, ref_fasta, gff3_file, vcf_file, tmp_path):
        """--mode pool + --samples → error."""
        result = runner.invoke(main, [
            "run", "--mode", "pool", "--samples", "S1,S2",
            "--vcf", vcf_file, "--gff", gff3_file,
            "--fasta", ref_fasta, "--outdir", str(tmp_path),
        ])
        assert result.exit_code != 0
        assert "individual" in result.output.lower() or "pool" in result.output.lower()

    def test_pool_with_samples_file_error(self, runner, ref_fasta, gff3_file, vcf_file, tmp_path):
        """--mode pool + --samples-file → error."""
        sf = tmp_path / "samples.txt"
        sf.write_text("S1\nS2\n")
        result = runner.invoke(main, [
            "run", "--mode", "pool", "--samples-file", str(sf),
            "--vcf", vcf_file, "--gff", gff3_file,
            "--fasta", ref_fasta, "--outdir", str(tmp_path),
        ])
        assert result.exit_code != 0

    def test_pool_with_min_call_rate_error(self, runner, ref_fasta, gff3_file, vcf_file, tmp_path):
        """--mode pool + --min-call-rate → error."""
        result = runner.invoke(main, [
            "run", "--mode", "pool", "--min-call-rate", "0.5",
            "--vcf", vcf_file, "--gff", gff3_file,
            "--fasta", ref_fasta, "--outdir", str(tmp_path),
        ])
        assert result.exit_code != 0

    def test_pool_with_min_an_error(self, runner, ref_fasta, gff3_file, vcf_file, tmp_path):
        """--mode pool + --min-an → error."""
        result = runner.invoke(main, [
            "run", "--mode", "pool", "--min-an", "4",
            "--vcf", vcf_file, "--gff", gff3_file,
            "--fasta", ref_fasta, "--outdir", str(tmp_path),
        ])
        assert result.exit_code != 0

    def test_individual_with_sample_error(self, runner, ref_fasta, gff3_file, vcf_file, tmp_path):
        """--mode individual + --sample → error."""
        result = runner.invoke(main, [
            "run", "--mode", "individual", "--sample", "S1",
            "--vcf", vcf_file, "--gff", gff3_file,
            "--fasta", ref_fasta, "--outdir", str(tmp_path),
        ])
        assert result.exit_code != 0

    def test_individual_with_min_depth_error(self, runner, ref_fasta, gff3_file, vcf_file, tmp_path):
        """--mode individual + --min-depth → error."""
        result = runner.invoke(main, [
            "run", "--mode", "individual", "--min-depth", "10",
            "--vcf", vcf_file, "--gff", gff3_file,
            "--fasta", ref_fasta, "--outdir", str(tmp_path),
        ])
        assert result.exit_code != 0

    def test_samples_and_samples_file_mutually_exclusive(self, runner, ref_fasta, gff3_file,
                                                          individual_vcf_file, tmp_path):
        sf = tmp_path / "samples.txt"
        sf.write_text("S1\nS2\n")
        result = runner.invoke(main, [
            "run", "--mode", "individual",
            "--samples", "S1,S2", "--samples-file", str(sf),
            "--vcf", individual_vcf_file, "--gff", gff3_file,
            "--fasta", ref_fasta, "--outdir", str(tmp_path),
        ])
        assert result.exit_code != 0
        assert "mutually exclusive" in result.output.lower()

    def test_min_call_rate_out_of_range(self, runner, ref_fasta, gff3_file, vcf_file, tmp_path):
        result = runner.invoke(main, [
            "run", "--mode", "individual", "--min-call-rate", "1.5",
            "--vcf", vcf_file, "--gff", gff3_file,
            "--fasta", ref_fasta, "--outdir", str(tmp_path),
        ])
        assert result.exit_code != 0

    def test_individual_no_samples_uses_all(self, runner, ref_fasta, gff3_file,
                                             individual_vcf_file, tmp_path):
        """--mode individual without --samples → uses all samples in VCF."""
        result = runner.invoke(main, [
            "run", "--mode", "individual",
            "--vcf", individual_vcf_file, "--gff", gff3_file,
            "--fasta", ref_fasta, "--outdir", str(tmp_path),
            "--min-freq", "0", "--min-qual", "0", "--min-call-rate", "0",
        ])
        assert result.exit_code == 0, result.output

    def test_invalid_sample_name_error(self, runner, ref_fasta, gff3_file,
                                       individual_vcf_file, tmp_path):
        """--samples with name not in VCF → error with available names."""
        result = runner.invoke(main, [
            "run", "--mode", "individual",
            "--samples", "S1,NONEXISTENT",
            "--vcf", individual_vcf_file, "--gff", gff3_file,
            "--fasta", ref_fasta, "--outdir", str(tmp_path),
        ])
        assert result.exit_code != 0
        assert "NONEXISTENT" in result.output

    def test_samples_file(self, runner, ref_fasta, gff3_file,
                          individual_vcf_file, tmp_path):
        """--samples-file reads sample names from file."""
        sf = tmp_path / "samples.txt"
        sf.write_text("S1\nS2\n")
        result = runner.invoke(main, [
            "run", "--mode", "individual",
            "--samples-file", str(sf),
            "--vcf", individual_vcf_file, "--gff", gff3_file,
            "--fasta", ref_fasta, "--outdir", str(tmp_path),
            "--min-freq", "0", "--min-qual", "0", "--min-call-rate", "0",
        ])
        assert result.exit_code == 0, result.output
```

**Step 2: Run tests to verify they fail**

Run: `mamba run -n pie pytest tests/test_cli.py::TestModeValidation -v`
Expected: FAIL (no `--mode` option yet).

**Step 3: Implement CLI changes**

Replace the `run` command in `src/pie/cli.py`:

```python
@main.command()
@click.option("--vcf", required=True, help="Input VCF file (bgzipped or plain)")
@click.option("--gff", required=True, help="GFF3 or GTF annotation file")
@click.option("--fasta", required=True, help="Reference FASTA file (indexed)")
@click.option("--outdir", required=True, help="Output directory")
@click.option("--mode", default="pool", show_default=True,
              help="Analysis mode: 'pool' for pool-seq (AD-based), "
                   "'individual' for individual sequencing (GT-based). "
                   "Alias 'ind' accepted for 'individual'.")
@click.option("--min-freq", default=0.01, show_default=True, help="Minimum allele frequency")
@click.option("--min-depth", default=None, type=int,
              help="Minimum read depth [pool mode only, default: 10]")
@click.option("--min-qual", default=20.0, show_default=True, help="Minimum variant quality")
@click.option("--pass-only", is_flag=True, help="Only use PASS-filtered variants")
@click.option("--keep-multiallelic", is_flag=True,
              help="Keep and merge multiallelic sites instead of skipping them")
@click.option("--include-stop-codons", is_flag=True,
              help="Count stop_gained mutations as nonsynonymous "
                   "(by default they are excluded, matching NG86/SNPGenie conventions)")
@click.option("--window-size", default=1000, show_default=True, help="Sliding window size (bp)")
@click.option("--window-step", default=100, show_default=True,
              type=click.IntRange(min=1), help="Sliding window step (bp)")
@click.option("--threads", default=1, show_default=True, help="Number of threads")
@click.option("--sample", default=None,
              help="Sample name to analyse [pool mode only]")
@click.option("--samples", default=None,
              help="Comma-separated sample names [individual mode only]")
@click.option("--samples-file", default=None, type=click.Path(exists=True),
              help="File with one sample name per line [individual mode only]")
@click.option("--min-call-rate", default=None, type=float,
              help="Minimum genotype call rate [individual mode only, default: 0.8]")
@click.option("--min-an", default=None, type=int,
              help="Minimum allele number (AN) [individual mode only, default: 2]")
def run(vcf, gff, fasta, outdir, mode, min_freq, min_depth, min_qual,
        pass_only, keep_multiallelic, include_stop_codons, window_size,
        window_step, threads, sample, samples, samples_file,
        min_call_rate, min_an):
    """Run piN/piS analysis."""
    from pie.vcf import ensure_indexed, get_sample_names
    from pie.parallel import run_parallel

    # Normalize mode alias
    if mode == "ind":
        mode = "individual"
    if mode not in ("pool", "individual"):
        click.echo(f"Error: --mode must be 'pool' or 'individual' (got '{mode}')", err=True)
        sys.exit(1)

    # --- Cross-option validation ---
    pool_only_given = {
        "--sample": sample is not None,
        "--min-depth": min_depth is not None,
    }
    ind_only_given = {
        "--samples": samples is not None,
        "--samples-file": samples_file is not None,
        "--min-call-rate": min_call_rate is not None,
        "--min-an": min_an is not None,
    }

    if mode == "pool":
        for opt, given in ind_only_given.items():
            if given:
                click.echo(
                    f"Error: {opt} is only valid in individual mode (--mode individual)",
                    err=True,
                )
                sys.exit(1)
        # Apply pool defaults
        if min_depth is None:
            min_depth = 10
    else:
        for opt, given in pool_only_given.items():
            if given:
                click.echo(
                    f"Error: {opt} is only valid in pool mode (--mode pool)",
                    err=True,
                )
                sys.exit(1)
        # Apply individual defaults
        if min_call_rate is None:
            min_call_rate = 0.8
        if min_an is None:
            min_an = 2
        if min_depth is None:
            min_depth = 10  # unused but needed for run_parallel signature

    # --samples and --samples-file mutually exclusive
    if samples is not None and samples_file is not None:
        click.echo(
            "Error: --samples and --samples-file are mutually exclusive",
            err=True,
        )
        sys.exit(1)

    # --min-call-rate range
    if min_call_rate is not None and not (0.0 <= min_call_rate <= 1.0):
        click.echo(
            f"Error: --min-call-rate must be between 0 and 1 (got {min_call_rate})",
            err=True,
        )
        sys.exit(1)

    # Validate inputs exist
    for path, name in [(vcf, "VCF"), (gff, "GFF"), (fasta, "FASTA")]:
        if not os.path.exists(path):
            click.echo(f"Error: {name} file not found: {path}", err=True)
            sys.exit(1)

    # Resolve sample list
    vcf_samples = get_sample_names(vcf)

    if mode == "pool":
        # Pool mode: existing sample validation
        if sample is not None:
            if sample not in vcf_samples:
                click.echo(
                    f"Error: sample '{sample}' not found in VCF. "
                    f"Available samples: {', '.join(vcf_samples)}",
                    err=True,
                )
                sys.exit(1)
            log.info("Using sample: %s", sample)
        elif len(vcf_samples) >= 2:
            click.echo(
                f"Error: VCF contains {len(vcf_samples)} samples. "
                f"Use --sample to specify one.\n"
                f"Available samples: {', '.join(vcf_samples)}",
                err=True,
            )
            sys.exit(1)
        selected_samples = None
    else:
        # Individual mode: resolve selected samples
        if samples is not None:
            selected_samples = [s.strip() for s in samples.split(",") if s.strip()]
        elif samples_file is not None:
            with open(samples_file) as f:
                selected_samples = [line.strip() for line in f if line.strip()]
        else:
            selected_samples = list(vcf_samples)

        # Validate sample names
        missing = [s for s in selected_samples if s not in vcf_samples]
        if missing:
            click.echo(
                f"Error: samples not found in VCF: {', '.join(missing)}\n"
                f"Available samples: {', '.join(vcf_samples)}",
                err=True,
            )
            sys.exit(1)
        log.info("Individual mode: %d samples selected", len(selected_samples))

    # Create output directory
    os.makedirs(outdir, exist_ok=True)

    # Ensure VCF is indexed
    vcf = ensure_indexed(vcf)
    log.info("VCF ready: %s", vcf)

    # Run analysis
    log.info("Starting piN/piS analysis (%s mode) with %d thread(s)", mode, threads)
    results = run_parallel(
        fasta_path=fasta, gff_path=gff, vcf_path=vcf,
        min_freq=min_freq, min_depth=min_depth, min_qual=min_qual,
        pass_only=pass_only, keep_multiallelic=keep_multiallelic,
        exclude_stops=not include_stop_codons, threads=threads,
        sample=sample, mode=mode, samples=selected_samples,
        min_call_rate=min_call_rate, min_an=min_an,
    )
    log.info("Processed %d genes", len(results))

    # Write outputs
    from pie.io import write_gene_results, write_window_results, write_summary
    prefix = f"{sample}." if sample else ""

    gene_path = os.path.join(outdir, f"{prefix}gene_results.tsv")
    write_gene_results(results, gene_path)
    log.info("Gene results: %s", gene_path)

    if window_size > 0:
        win_path = os.path.join(outdir, f"{prefix}window_results.tsv")
        write_window_results(results, win_path, window_size, window_step)
        log.info("Window results: %s", win_path)

    summary_path = os.path.join(outdir, f"{prefix}summary.tsv")
    write_summary(results, summary_path)
    log.info("Summary: %s", summary_path)

    click.echo(f"Done. Results written to {outdir}/")
```

**Step 4: Run CLI validation tests**

Run: `mamba run -n pie pytest tests/test_cli.py -v`
Expected: All PASS.

**Step 5: Run existing tests for pool-mode regression**

Run: `mamba run -n pie pytest tests/ -v`
Expected: All PASS.

**Step 6: Commit**

```bash
git add src/pie/cli.py tests/test_cli.py
git commit -m "feat: add --mode pool|individual CLI with cross-option validation"
```

---

### Task 10: End-to-end integration test for individual mode

**Files:**
- Modify: `tests/test_integration.py`

**Step 1: Write integration test with hand-calculated values**

Hand-calculated expected values for `individual_vcf_file` (all 4 samples, all filters disabled):

- **Gene1** (pos 1–90, + strand, 29 coding codons excl. stop):
  - pos 6 T→C (GCT→GCC, syn), freq = 2/6 = 1/3
  - pos 7 G→A (GAT→AAT, nonsyn), freq = 2/8 = 1/4
  - S_diffs(codon 2) = 2 × (2/3) × (1/3) = 4/9 ≈ 0.4444
  - N_diffs(codon 3) = 2 × (3/4) × (1/4) = 3/8 = 0.375
  - N_sites = 59.6667, S_sites = 27.3333 (same codons as pool test)

- **Gene2** (pos 101–220, + strand, 32 codons):
  - pos 195 A→T (GAT→GTT, nonsyn), freq = 6/8 = 3/4
  - N_diffs = 2 × (1/4) × (3/4) = 3/8 = 0.375
  - S_diffs = 0

- **Gene3** (pos 231–311, − strand):
  - pos 297 A→G (fwd), T→C sense, GAT→GAC (syn), freq = 1/4
  - S_diffs = 2 × (3/4) × (1/4) = 3/8 = 0.375
  - N_diffs = 0

Add to `tests/test_integration.py`:

```python
class TestIndividualMode:
    """End-to-end tests for --mode individual using GT-derived frequencies."""

    def _run_individual(self, runner, ref_fasta, gff3_file, individual_vcf_file,
                        tmp_path, min_qual=0, min_call_rate=0):
        result = runner.invoke(main, [
            "run", "--mode", "individual",
            "--vcf", individual_vcf_file,
            "--gff", gff3_file,
            "--fasta", ref_fasta,
            "--outdir", str(tmp_path),
            "--min-freq", "0",
            "--min-qual", str(min_qual),
            "--min-call-rate", str(min_call_rate),
            "--min-an", "0",
            "--window-size", "30",
            "--window-step", "10",
        ])
        assert result.exit_code == 0, result.output
        return tmp_path

    def test_gene1_individual(self, runner, ref_fasta, gff3_file,
                               individual_vcf_file, tmp_path):
        """Gene1: pos6 T>C syn freq=1/3, pos7 G>A nonsyn freq=1/4."""
        out = self._run_individual(runner, ref_fasta, gff3_file,
                                    individual_vcf_file, tmp_path)
        df = pd.read_csv(out / "gene_results.tsv", sep="\t")
        g1 = df[df["gene_id"].str.contains("gene1", case=False)].iloc[0]

        assert g1["n_codons"] == 29
        assert g1["n_poly_codons"] == 2
        assert g1["n_variants"] == 2
        assert abs(g1["N_sites"] - 59.6667) < 0.001
        assert abs(g1["S_sites"] - 27.3333) < 0.001
        assert abs(g1["N_diffs"] - 0.375) < 1e-6    # 2*(3/4)*(1/4)
        assert abs(g1["S_diffs"] - 4 / 9) < 1e-6    # 2*(2/3)*(1/3)

    def test_gene2_individual(self, runner, ref_fasta, gff3_file,
                               individual_vcf_file, tmp_path):
        """Gene2: pos195 A>T nonsyn freq=3/4."""
        out = self._run_individual(runner, ref_fasta, gff3_file,
                                    individual_vcf_file, tmp_path)
        df = pd.read_csv(out / "gene_results.tsv", sep="\t")
        g2 = df[df["gene_id"].str.contains("gene2", case=False)].iloc[0]

        assert g2["n_poly_codons"] == 1
        assert abs(g2["N_diffs"] - 0.375) < 1e-6    # 2*(1/4)*(3/4)
        assert g2["S_diffs"] == 0.0

    def test_gene3_individual(self, runner, ref_fasta, gff3_file,
                               individual_vcf_file, tmp_path):
        """Gene3: pos297 syn freq=1/4."""
        out = self._run_individual(runner, ref_fasta, gff3_file,
                                    individual_vcf_file, tmp_path)
        df = pd.read_csv(out / "gene_results.tsv", sep="\t")
        g3 = df[df["gene_id"].str.contains("gene3", case=False)].iloc[0]

        assert abs(g3["S_diffs"] - 0.375) < 1e-6    # 2*(3/4)*(1/4)
        assert g3["N_diffs"] == 0.0

    def test_output_has_sample_metadata(self, runner, ref_fasta, gff3_file,
                                         individual_vcf_file, tmp_path):
        """Individual mode outputs include n_samples and mean_call_rate."""
        out = self._run_individual(runner, ref_fasta, gff3_file,
                                    individual_vcf_file, tmp_path)
        df = pd.read_csv(out / "gene_results.tsv", sep="\t")
        assert "n_samples" in df.columns
        assert "mean_call_rate" in df.columns
        assert (df["n_samples"] == 4).all()

    def test_summary_has_sample_metadata(self, runner, ref_fasta, gff3_file,
                                          individual_vcf_file, tmp_path):
        out = self._run_individual(runner, ref_fasta, gff3_file,
                                    individual_vcf_file, tmp_path)
        df = pd.read_csv(out / "summary.tsv", sep="\t")
        assert "n_samples_selected" in df.columns
        assert df.iloc[0]["n_samples_selected"] == 4

    def test_min_call_rate_filters_variants(self, runner, ref_fasta, gff3_file,
                                             individual_vcf_file, tmp_path):
        """min_call_rate=0.8 filters pos6 (0.75) and pos297 (0.50)."""
        out = self._run_individual(runner, ref_fasta, gff3_file,
                                    individual_vcf_file, tmp_path,
                                    min_call_rate=0.8)
        df = pd.read_csv(out / "gene_results.tsv", sep="\t")
        g1 = df[df["gene_id"].str.contains("gene1", case=False)].iloc[0]
        # Only pos7 passes → 1 variant in gene1
        assert g1["n_variants"] == 1

    def test_all_output_files_created(self, runner, ref_fasta, gff3_file,
                                       individual_vcf_file, tmp_path):
        out = self._run_individual(runner, ref_fasta, gff3_file,
                                    individual_vcf_file, tmp_path)
        assert (out / "gene_results.tsv").exists()
        assert (out / "window_results.tsv").exists()
        assert (out / "summary.tsv").exists()
```

**Step 2: Run integration tests**

Run: `mamba run -n pie pytest tests/test_integration.py::TestIndividualMode -v`
Expected: All PASS.

**Step 3: Commit**

```bash
git add tests/test_integration.py
git commit -m "test: add end-to-end integration tests for individual mode"
```

---

### Task 11: Pool mode regression — full test suite

**Files:** none (verification only)

**Step 1: Run the complete test suite**

Run: `mamba run -n pie pytest tests/ -v`
Expected: All existing tests PASS. No pool-mode behavior changed.

**Step 2: Verify pool-mode output format unchanged**

Run: `mamba run -n pie pytest tests/test_integration.py::TestEndToEnd -v`
Expected: All PASS. No `n_samples` or `mean_call_rate` columns in pool-mode output.

**Step 3: Commit (only if any fixup needed)**

If all green, no commit needed. Tag as final verification.

---

### Task 12: Update README and docstring

**Files:**
- Modify: `src/pie/cli.py:1` (module docstring)

**Step 1: Update module docstring**

```python
"""pie CLI — piN/piS Estimator for pool-seq and individual-sequencing data."""
```

**Step 2: Commit**

```bash
git add src/pie/cli.py
git commit -m "docs: update CLI docstring for individual mode support"
```
