"""Core piN/piS diversity engine using the Nei-Gojobori method.

Computes per-gene nonsynonymous (piN) and synonymous (piS) nucleotide
diversity from pooled sequencing data.
"""

import logging
from dataclasses import dataclass, field
from itertools import product

import numpy as np

from pie.codon import (
    CODON_TO_INDEX,
    N_DIFFS,
    N_SITES,
    S_DIFFS,
    S_SITES,
    AMINO_ACID,
    is_stop_codon,
)
from pie.annotation import GeneModel
from pie.reference import ReferenceGenome
from pie.vcf import VariantReader, Variant

log = logging.getLogger(__name__)

# Base encoding: A=0, C=1, G=2, T=3
_BASE_TO_IDX = {"A": 0, "C": 1, "G": 2, "T": 3}
_IDX_TO_BASE = "ACGT"
_COMPLEMENT_BASE = {"A": "T", "T": "A", "C": "G", "G": "C"}


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
            freq[i, j, _BASE_TO_IDX[base]] = 1.0

    if not variants:
        return freq

    # Build lookup: genomic_pos -> (codon_idx, pos_within_codon)
    pos_map: dict[int, tuple[int, int]] = {}
    for i, (chrom, p1, p2, p3) in enumerate(positions):
        pos_map[p1] = (i, 0)
        pos_map[p2] = (i, 1)
        pos_map[p3] = (i, 2)

    # Overlay variants
    for var in variants:
        if var.pos not in pos_map:
            continue
        codon_idx, pos_in_codon = pos_map[var.pos]

        # Convert variant alleles to CDS-sense
        alt = var.alt
        if strand == "-":
            alt = _COMPLEMENT_BASE[alt]

        alt_idx = _BASE_TO_IDX[alt]

        # Reduce reference frequency, add alt frequency
        # The reference allele in the freq array is already CDS-sense
        ref_base = codons[codon_idx][pos_in_codon]
        ref_idx = _BASE_TO_IDX[ref_base]

        freq[codon_idx, pos_in_codon, ref_idx] -= var.freq
        freq[codon_idx, pos_in_codon, alt_idx] += var.freq

    # Clamp to [0, 1] and renormalize
    freq = np.clip(freq, 0.0, 1.0)
    sums = freq.sum(axis=2, keepdims=True)
    sums[sums == 0] = 1.0  # avoid division by zero
    freq /= sums

    return freq


# ---------------------------------------------------------------------------
# 2. Compute diversity for a single codon
# ---------------------------------------------------------------------------
def compute_codon_diversity(freq: np.ndarray) -> dict:
    """Compute N/S sites and diffs for a single codon.

    Args:
        freq: shape (3, 4) allele frequency array for one codon.

    Returns:
        dict with keys N_sites, S_sites, N_diffs, S_diffs.
    """
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

    # Handle stop codons: remove and renormalize
    stop_freq = sum(p for _, p, idx in codon_probs if is_stop_codon(idx))
    if stop_freq > 0.01:
        log.warning(
            "Stop codon frequency %.4f > 1%% in polymorphic codon; "
            "removing and renormalizing",
            stop_freq,
        )

    if stop_freq > 0:
        codon_probs = [(s, p, i) for s, p, i in codon_probs if not is_stop_codon(i)]
        total = sum(p for _, p, _ in codon_probs)
        if total > 0:
            codon_probs = [(s, p / total, i) for s, p, i in codon_probs]
        else:
            # All codons are stops (shouldn't happen); return zeros
            return {"N_sites": 0.0, "S_sites": 0.0, "N_diffs": 0.0, "S_diffs": 0.0}

    # Weighted site counts: N_sites = sum(freq_i * N_SITES[codon_i].sum())
    n_sites = 0.0
    s_sites = 0.0
    for _, prob, idx in codon_probs:
        n_sites += prob * float(N_SITES[idx].sum())
        s_sites += prob * float(S_SITES[idx].sum())

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
    }


# ---------------------------------------------------------------------------
# 3. Compute diversity for a whole gene
# ---------------------------------------------------------------------------
def compute_gene_diversity(
    gene: GeneModel,
    ref: ReferenceGenome,
    vcf: VariantReader,
) -> GeneResult:
    """Compute per-gene piN/piS diversity.

    Args:
        gene: GeneModel with CDS exon coordinates.
        ref: ReferenceGenome for sequence access.
        vcf: VariantReader for variant access.

    Returns:
        GeneResult with accumulated N/S sites and diffs.
    """
    # Extract codons and genomic positions
    codons = ref.extract_codons(gene.cds_exons, gene.strand)
    positions = ref.codon_genomic_positions(gene.cds_exons, gene.strand)

    # Fetch all variants across CDS exons
    all_variants: list[Variant] = []
    for chrom, start, end in gene.cds_exons:
        all_variants.extend(vcf.fetch(chrom, start, end))

    # Build allele frequency array
    freq_array = build_allele_freq_array(codons, positions, all_variants, gene.strand)

    # Identify polymorphic codons: any position with >1 allele having freq > 0
    poly_mask = np.zeros(len(codons), dtype=bool)
    for i in range(len(codons)):
        for j in range(3):
            n_alleles = np.sum(freq_array[i, j] > 0)
            if n_alleles > 1:
                poly_mask[i] = True
                break

    # Accumulate results
    total_N_sites = 0.0
    total_S_sites = 0.0
    total_N_diffs = 0.0
    total_S_diffs = 0.0
    n_codons_analyzed = 0
    n_poly = 0
    codon_results: list[CodonResult] = []

    for i, codon_str in enumerate(codons):
        # Skip stop codons entirely
        idx = CODON_TO_INDEX.get(codon_str)
        if idx is not None and is_stop_codon(idx):
            continue

        chrom = positions[i][0]
        pos1 = positions[i][1]
        n_codons_analyzed += 1

        if poly_mask[i]:
            # Full diversity computation for polymorphic codons
            n_poly += 1
            result = compute_codon_diversity(freq_array[i])
            cr = CodonResult(
                chrom=chrom,
                pos1=pos1,
                N_sites=result["N_sites"],
                S_sites=result["S_sites"],
                N_diffs=result["N_diffs"],
                S_diffs=result["S_diffs"],
            )
        else:
            # Monomorphic: only site counts, diffs = 0
            if idx is not None:
                n_s = float(N_SITES[idx].sum())
                s_s = float(S_SITES[idx].sum())
            else:
                n_s = 0.0
                s_s = 0.0
            cr = CodonResult(
                chrom=chrom,
                pos1=pos1,
                N_sites=n_s,
                S_sites=s_s,
                N_diffs=0.0,
                S_diffs=0.0,
            )

        codon_results.append(cr)
        total_N_sites += cr.N_sites
        total_S_sites += cr.S_sites
        total_N_diffs += cr.N_diffs
        total_S_diffs += cr.S_diffs

    # Mean depth from variants
    if all_variants:
        mean_depth = sum(v.depth for v in all_variants) / len(all_variants)
    else:
        mean_depth = 0.0

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
        mean_depth=mean_depth,
        n_variants=len(all_variants),
        codon_results=codon_results,
    )
