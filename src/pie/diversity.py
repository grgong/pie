"""Core piN/piS diversity engine using the Nei-Gojobori method.

Computes per-gene nonsynonymous (piN) and synonymous (piS) nucleotide
diversity from pool-seq or individual-sequencing data.
"""

import logging
from dataclasses import dataclass, field
from itertools import product

import numpy as np

from pie.codon import (
    CODON_TO_INDEX,
    N_DIFFS,
    N_SITES,
    N_SITES_EXCL_STOP,
    N_SITES_EXCL_STOP_SUM,
    N_SITES_SUM,
    S_DIFFS,
    S_SITES,
    S_SITES_EXCL_STOP,
    S_SITES_EXCL_STOP_SUM,
    S_SITES_SUM,
    is_stop_codon,
)
from typing import Protocol

from pie.annotation import GeneModel
from pie.reference import ReferenceGenome
from pie.vcf import FetchResult, FilterStats, Variant

log = logging.getLogger(__name__)


class VariantReaderLike(Protocol):
    """Protocol for variant readers (VariantReader or IndividualVariantReader)."""

    def fetch(self, chrom: str, start: int, end: int) -> FetchResult: ...

# Base encoding: A=0, C=1, G=2, T=3
_BASE_TO_IDX = {"A": 0, "C": 1, "G": 2, "T": 3}
_IDX_TO_BASE = "ACGT"
_COMPLEMENT_BASE = {"A": "T", "T": "A", "C": "G", "G": "C"}
_VALID_CODON_BASES = frozenset("ACGT")


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------
@dataclass
class CodonResult:
    """Per-codon diversity result."""

    chrom: str
    pos1: int  # genomic position of first codon base (0-based)
    N_sites: float
    S_sites: float
    N_diffs: float
    S_diffs: float


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
    n_stop_codons: int = 0
    codon_results: list[CodonResult] = field(default_factory=list)
    n_samples: int | None = None
    call_rates: list[float] | None = None
    # QC / filter statistics
    n_ambiguous_codons: int = 0
    n_internal_stop_codons: int = 0
    filter_stats: FilterStats = field(default_factory=FilterStats)

    @property
    def mean_call_rate(self) -> float | None:
        if self.call_rates is None or len(self.call_rates) == 0:
            return None
        return sum(self.call_rates) / len(self.call_rates)

    @property
    def piN(self) -> float:
        return self.N_diffs / self.N_sites if self.N_sites > 0 else 0.0

    @property
    def piS(self) -> float:
        return self.S_diffs / self.S_sites if self.S_sites > 0 else 0.0

    @property
    def piN_piS(self) -> float | None:
        return self.piN / self.piS if self.piS > 0 else None


# ---------------------------------------------------------------------------
# 1. Build allele frequency array
# ---------------------------------------------------------------------------
def build_allele_freq_array(
    codons: list[str],
    positions: list[tuple[str, int, int, int]],
    variants: list[Variant],
    strand: str = "+",
) -> np.ndarray:
    """Build a (n_codons, 3, 4) array of allele frequencies.

    Args:
        codons: CDS-sense codon strings (already reverse-complemented for - strand).
        positions: Genomic positions per codon [(chrom, pos1, pos2, pos3), ...].
        variants: List of Variant objects with forward-strand coordinates.
        strand: "+" or "-". Determines whether variant alleles need complementing.

    Returns:
        np.ndarray of shape (n_codons, 3, 4) with allele frequencies.
    """
    n = len(codons)
    freq = np.zeros((n, 3, 4), dtype=np.float64)

    # Initialize with reference alleles at frequency 1.0
    for i, codon in enumerate(codons):
        for j, base in enumerate(codon):
            idx = _BASE_TO_IDX.get(base)
            if idx is None:
                # Non-ACGT base (e.g. N) — leave row as zeros; upstream
                # should have filtered this codon, but guard defensively.
                break
            freq[i, j, idx] = 1.0

    if not variants:
        return freq

    # Build lookup: genomic_pos -> (codon_idx, pos_within_codon)
    pos_map: dict[int, tuple[int, int]] = {}
    for i, (chrom, p1, p2, p3) in enumerate(positions):
        pos_map[p1] = (i, 0)
        pos_map[p2] = (i, 1)
        pos_map[p3] = (i, 2)

    # Group variants by position for proper multi-allelic handling
    pos_variants: dict[int, list[Variant]] = {}
    for var in variants:
        if var.pos in pos_map:
            pos_variants.setdefault(var.pos, []).append(var)

    # Apply variant frequencies per position (handles multi-allelic correctly)
    for pos, var_list in pos_variants.items():
        codon_idx, pos_in_codon = pos_map[pos]
        ref_base = codons[codon_idx][pos_in_codon]
        ref_idx = _BASE_TO_IDX[ref_base]

        # Set each alt allele frequency directly
        total_alt_freq = 0.0
        for var in var_list:
            alt = var.alt
            if strand == "-":
                alt = _COMPLEMENT_BASE[alt]
            alt_idx = _BASE_TO_IDX[alt]
            freq[codon_idx, pos_in_codon, alt_idx] = var.freq
            total_alt_freq += var.freq

        # Reference = 1 - sum(alt frequencies)
        freq[codon_idx, pos_in_codon, ref_idx] = max(0.0, 1.0 - total_alt_freq)

    # Clamp to [0, 1] and renormalize
    freq = np.clip(freq, 0.0, 1.0)
    sums = freq.sum(axis=2, keepdims=True)
    sums[sums == 0] = 1.0  # avoid division by zero
    freq /= sums

    return freq


# ---------------------------------------------------------------------------
# 2. Compute diversity for a single codon
# ---------------------------------------------------------------------------
def compute_codon_diversity(freq: np.ndarray, exclude_stops: bool = False) -> dict:
    """Compute N/S sites and diffs for a single codon.

    Args:
        freq: shape (3, 4) allele frequency array for one codon.
        exclude_stops: If True (default at CLI level), remove stop codon
            alleles and renormalize (NG86/SNPGenie convention). If False,
            keep stop alleles and count them as nonsynonymous.

    Returns:
        dict with keys N_sites, S_sites, N_diffs, S_diffs.
    """
    # Choose site tables
    if exclude_stops:
        n_sites_tbl, s_sites_tbl = N_SITES_EXCL_STOP, S_SITES_EXCL_STOP
    else:
        n_sites_tbl, s_sites_tbl = N_SITES, S_SITES

    # Enumerate possible codons: alleles with freq > 0 at each position
    alleles_per_pos = []
    for pos in range(3):
        alleles = [(base_idx, freq[pos, base_idx])
                   for base_idx in range(4)
                   if freq[pos, base_idx] > 0]
        alleles_per_pos.append(alleles)

    # Build all possible codons with their probabilities
    codon_probs: list[tuple[str, float, int]] = []  # (codon_str, prob, codon_index)
    for combo in product(*alleles_per_pos):
        bases = [_IDX_TO_BASE[c[0]] for c in combo]
        prob = combo[0][1] * combo[1][1] * combo[2][1]
        codon_str = "".join(bases)
        codon_idx = CODON_TO_INDEX[codon_str]
        codon_probs.append((codon_str, prob, codon_idx))

    # Handle stop codons
    stop_freq = 0.0
    if exclude_stops:
        stop_freq = sum(p for _, p, idx in codon_probs if is_stop_codon(idx))

        if stop_freq > 0:
            codon_probs = [(s, p, i) for s, p, i in codon_probs if not is_stop_codon(i)]
            total = sum(p for _, p, _ in codon_probs)
            if total > 0:
                codon_probs = [(s, p / total, i) for s, p, i in codon_probs]
            else:
                return {"N_sites": 0.0, "S_sites": 0.0, "N_diffs": 0.0, "S_diffs": 0.0,
                        "_stop_freq": stop_freq}

    # Weighted site counts
    n_sites = 0.0
    s_sites = 0.0
    for _, prob, idx in codon_probs:
        n_sites += prob * float(n_sites_tbl[idx].sum())
        s_sites += prob * float(s_sites_tbl[idx].sum())

    # Weighted pairwise diffs
    n_diffs = 0.0
    s_diffs = 0.0
    n_codons = len(codon_probs)
    for i in range(n_codons):
        for j in range(i + 1, n_codons):
            _, pi, idx_i = codon_probs[i]
            _, pj, idx_j = codon_probs[j]
            weight = 2.0 * pi * pj
            n_diffs += weight * float(N_DIFFS[idx_i, idx_j])
            s_diffs += weight * float(S_DIFFS[idx_i, idx_j])

    return {
        "N_sites": n_sites,
        "S_sites": s_sites,
        "N_diffs": n_diffs,
        "S_diffs": s_diffs,
        "_stop_freq": stop_freq,
    }


def _monomorphic_codon_index(freq: np.ndarray) -> int | None:
    """Return the CODON_TO_INDEX for the actual codon encoded in *freq*.

    For a monomorphic codon each position has exactly one allele with
    frequency > 0.  This function reads those alleles and returns the
    codon index, which may differ from the reference when the population
    is fixed for an alternate allele.

    Returns None if any position has no allele with freq > 0.
    """
    bases = []
    for pos in range(3):
        nonzero = np.nonzero(freq[pos])[0]
        if len(nonzero) != 1:
            return None
        bases.append(_IDX_TO_BASE[nonzero[0]])
    return CODON_TO_INDEX.get("".join(bases))


# ---------------------------------------------------------------------------
# 3. Compute diversity for a whole gene
# ---------------------------------------------------------------------------
def compute_gene_diversity(
    gene: GeneModel,
    ref: ReferenceGenome,
    vcf: VariantReaderLike,
    exclude_stops: bool = False,
) -> GeneResult:
    """Compute per-gene piN/piS diversity.

    Args:
        gene: GeneModel with CDS exon coordinates.
        ref: ReferenceGenome for sequence access.
        vcf: VariantReader for variant access.
        exclude_stops: If True (default at CLI level), exclude stop_gained
            mutations from piN (NG86/SNPGenie convention). If False,
            count them as nonsynonymous.

    Returns:
        GeneResult with accumulated N/S sites and diffs.
    """
    # Extract codons and genomic positions
    codons = ref.extract_codons(gene.cds_exons, gene.strand)
    positions = ref.codon_genomic_positions(gene.cds_exons, gene.strand)

    # Filter out codons containing non-ACGT bases (e.g. N in reference)
    n_ambiguous = 0
    if any(base not in _VALID_CODON_BASES for codon in codons for base in codon):
        clean = [(c, p) for c, p in zip(codons, positions)
                 if all(b in _VALID_CODON_BASES for b in c)]
        n_ambiguous = len(codons) - len(clean)
        if clean:
            codons, positions = zip(*clean)
            codons = list(codons)
            positions = list(positions)
        else:
            codons, positions = [], []
        log.debug("Gene %s: skipped %d codon(s) with ambiguous bases",
                  gene.gene_id, n_ambiguous)

    # Early return when no valid codons remain (all bases ambiguous)
    if not codons:
        return GeneResult(
            gene_id=gene.gene_id, transcript_id=gene.transcript_id,
            chrom=gene.chrom, start=gene.start, end=gene.end,
            strand=gene.strand, n_codons=0, n_poly_codons=0,
            N_sites=0.0, S_sites=0.0, N_diffs=0.0, S_diffs=0.0,
            mean_variant_depth=0.0, n_variants=0,
            n_ambiguous_codons=n_ambiguous,
        )

    # Fetch all variants across CDS exons
    all_variants: list[Variant] = []
    gene_filter_stats = FilterStats()
    for chrom, start, end in gene.cds_exons:
        result = vcf.fetch(chrom, start, end)
        all_variants.extend(result.variants)
        gene_filter_stats += result.stats

    # Choose site tables and precomputed row sums
    if exclude_stops:
        n_sites_tbl, s_sites_tbl = N_SITES_EXCL_STOP, S_SITES_EXCL_STOP
        n_sites_sum, s_sites_sum = N_SITES_EXCL_STOP_SUM, S_SITES_EXCL_STOP_SUM
    else:
        n_sites_tbl, s_sites_tbl = N_SITES, S_SITES
        n_sites_sum, s_sites_sum = N_SITES_SUM, S_SITES_SUM

    # --- Identify codons hit by at least one variant ---
    # Only these need the full freq_array + diversity computation;
    # the remaining ~98% use reference codon site counts directly.
    variant_positions = {v.pos for v in all_variants}
    hit_codon_map: dict[int, int] = {}
    hit_codons: list[str] = []
    hit_positions: list[tuple[str, int, int, int]] = []
    for i, (_, p1, p2, p3) in enumerate(positions):
        if p1 in variant_positions or p2 in variant_positions or p3 in variant_positions:
            hit_codon_map[i] = len(hit_codons)
            hit_codons.append(codons[i])
            hit_positions.append(positions[i])

    # Build allele frequency array only for variant-hit codons
    freq_array = None
    poly_mask = None
    if hit_codons:
        freq_array = build_allele_freq_array(
            hit_codons, hit_positions, all_variants, gene.strand)
        poly_mask = (freq_array > 0).sum(axis=2).max(axis=1) > 1

    # Accumulate results
    total_N_sites = 0.0
    total_S_sites = 0.0
    total_N_diffs = 0.0
    total_S_diffs = 0.0
    n_codons_analyzed = 0
    n_poly = 0
    n_stop_warn = 0  # count codons with stop freq > 1%
    n_internal_stops = 0  # count internal (non-terminal) stop codons
    codon_results: list[CodonResult] = []

    for i, codon_str in enumerate(codons):
        # Skip stop codons (always, regardless of mode)
        idx = CODON_TO_INDEX.get(codon_str)
        if idx is not None and is_stop_codon(idx):
            if i < len(codons) - 1:
                n_internal_stops += 1
            continue

        chrom = positions[i][0]
        pos1 = positions[i][1]
        n_codons_analyzed += 1

        j = hit_codon_map.get(i)
        if j is not None and poly_mask[j]:
            # Full diversity computation for polymorphic codons
            n_poly += 1
            div = compute_codon_diversity(freq_array[j], exclude_stops=exclude_stops)
            if exclude_stops and div.get("_stop_freq", 0.0) > 0.01:
                n_stop_warn += 1
            cr = CodonResult(
                chrom=chrom, pos1=pos1,
                N_sites=div["N_sites"], S_sites=div["S_sites"],
                N_diffs=div["N_diffs"], S_diffs=div["S_diffs"],
            )
        else:
            # Monomorphic codon: only site counts, diffs = 0.
            if j is not None:
                # Variant-hit but monomorphic (population fixed for alt allele)
                actual_idx = _monomorphic_codon_index(freq_array[j])
                if actual_idx is None:
                    actual_idx = idx
            else:
                # No variant touches this codon — use reference directly
                actual_idx = idx
            if actual_idx is not None:
                n_s = float(n_sites_sum[actual_idx])
                s_s = float(s_sites_sum[actual_idx])
            else:
                n_s = 0.0
                s_s = 0.0
            cr = CodonResult(
                chrom=chrom, pos1=pos1,
                N_sites=n_s, S_sites=s_s,
                N_diffs=0.0, S_diffs=0.0,
            )

        codon_results.append(cr)
        total_N_sites += cr.N_sites
        total_S_sites += cr.S_sites
        total_N_diffs += cr.N_diffs
        total_S_diffs += cr.S_diffs

    if n_internal_stops > 0:
        log.warning(
            "Gene %s: %d internal stop codon(s) in reference — "
            "possible pseudogene or annotation issue",
            gene.gene_id, n_internal_stops,
        )
    if n_stop_warn > 0:
        log.debug(
            "Gene %s: %d polymorphic codon(s) had stop-codon frequency > 1%%; "
            "stops removed and renormalized",
            gene.gene_id, n_stop_warn,
        )

    # Mean depth across variant sites only (0 when no variants)
    if all_variants:
        mean_variant_depth = sum(v.depth for v in all_variants) / len(all_variants)
    else:
        mean_variant_depth = 0.0

    # Collect per-variant call rates (individual mode only)
    cr_list = [v.call_rate for v in all_variants if v.call_rate is not None]

    return GeneResult(
        gene_id=gene.gene_id,
        transcript_id=gene.transcript_id,
        chrom=gene.chrom,
        start=gene.start,
        end=gene.end,
        strand=gene.strand,
        n_codons=n_codons_analyzed,
        n_poly_codons=n_poly,
        N_sites=total_N_sites,
        S_sites=total_S_sites,
        N_diffs=total_N_diffs,
        S_diffs=total_S_diffs,
        mean_variant_depth=mean_variant_depth,
        n_variants=len(all_variants),
        n_stop_codons=n_stop_warn,
        codon_results=codon_results,
        call_rates=cr_list if cr_list else None,
        n_ambiguous_codons=n_ambiguous,
        n_internal_stop_codons=n_internal_stops,
        filter_stats=gene_filter_stats,
    )
