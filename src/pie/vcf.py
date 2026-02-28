"""VCF parsing, filtering, and auto-indexing via cyvcf2."""
import os
import logging
from dataclasses import dataclass
from cyvcf2 import VCF
import pysam

log = logging.getLogger(__name__)

_VALID_BASES = frozenset("ACGT")


class _BaseVariantReader:
    """Shared context-manager, contig-tracking, and star-allele logging."""

    _vcf: VCF

    def _init_contig_tracking(self):
        """Call at end of subclass __init__ after self._vcf is set."""
        self._contigs: frozenset[str] = frozenset(self._vcf.seqnames)
        self._missing_contigs: set[str] = set()
        self._n_star_alleles = 0

    def close(self):
        if self._missing_contigs:
            log.warning(
                "Skipped %d contig(s) absent from VCF: %s",
                len(self._missing_contigs),
                ", ".join(sorted(self._missing_contigs)[:20])
                + (" ..." if len(self._missing_contigs) > 20 else ""),
            )
        if self._n_star_alleles:
            log.info("Skipped %d non-SNP ALT alleles (*, indels)",
                     self._n_star_alleles)
        self._vcf.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    def _check_contig(self, chrom: str) -> bool:
        """Return True if chrom exists in VCF; track missing ones."""
        if chrom not in self._contigs:
            self._missing_contigs.add(chrom)
            return False
        return True


@dataclass(slots=True)
class Variant:
    pos: int      # 0-based genomic position
    ref: str      # reference allele
    alt: str      # alternate allele
    freq: float   # alt allele frequency (0-1)
    depth: int    # total depth (read depth in pool mode, AN in individual mode)
    call_rate: float | None = None  # fraction of called samples (individual mode only)


def ensure_indexed(vcf_path: str) -> str:
    """Ensure VCF is bgzipped and tabix-indexed. Returns path to .vcf.gz."""
    need_index = False
    if vcf_path.endswith(".vcf"):
        gz_path = vcf_path + ".gz"
        log.info("Bgzipping %s -> %s", vcf_path, gz_path)
        pysam.tabix_compress(vcf_path, gz_path, force=True)
        vcf_path = gz_path
        need_index = True  # always re-index after fresh compression

    tbi_path = vcf_path + ".tbi"
    if need_index or not os.path.exists(tbi_path):
        log.info("Creating tabix index for %s", vcf_path)
        pysam.tabix_index(vcf_path, preset="vcf", force=True)

    return vcf_path


def get_sample_names(vcf_path: str) -> list[str]:
    """Return the list of sample names from a VCF file."""
    vcf = VCF(vcf_path)
    samples = list(vcf.samples)
    vcf.close()
    return samples


def get_vcf_contigs(vcf_path: str) -> frozenset[str]:
    """Return the set of contig/sequence names present in a VCF file."""
    vcf = VCF(vcf_path)
    contigs = frozenset(vcf.seqnames)
    vcf.close()
    return contigs


class VariantReader(_BaseVariantReader):
    def __init__(self, vcf_path: str, min_freq: float = 0.01,
                 min_depth: int = 10, min_qual: float = 20.0,
                 pass_only: bool = False, keep_multiallelic: bool = False,
                 sample: str | None = None):
        self._min_freq = min_freq
        self._min_depth = min_depth
        self._min_qual = min_qual
        self._pass_only = pass_only
        self._keep_multiallelic = keep_multiallelic
        if sample is not None:
            self._vcf = VCF(vcf_path, samples=[sample])
        else:
            self._vcf = VCF(vcf_path)
        self._init_contig_tracking()

    def fetch(self, chrom: str, start: int, end: int) -> list[Variant]:
        """Fetch filtered variants in region (0-based, half-open).

        Handles both single-line multiallelic records (``ALT=C,G``) and
        decomposed records (from ``bcftools norm -m-``).  By default,
        positions with multiple ALT alleles are skipped; set
        ``keep_multiallelic=True`` to merge them instead (frequencies
        recomputed from raw allele depths).
        """
        if not self._check_contig(chrom):
            return []

        region = f"{chrom}:{start + 1}-{end}"

        # Collect per-record data including raw allele depths
        # (pos, ref, alt, freq, depth, ref_count, alt_count)
        raw: list[tuple[int, str, str, float, int, int, int]] = []
        # Track positions with multiple ALT alleles BEFORE per-allele
        # filtering, so multiallelic detection is not affected by thresholds.
        multiallelic_pos: set[int] = set()
        seen_pos: dict[int, int] = {}  # pos -> count of valid SNP alleles
        for record in self._vcf(region):
            if self._pass_only and record.FILTER is not None:
                continue
            if record.QUAL is not None and record.QUAL < self._min_qual:
                continue
            # REF must be a single valid nucleotide for a SNP
            if record.REF not in _VALID_BASES:
                continue

            pos0 = record.POS - 1
            # Count valid SNP alleles at this position (before freq/depth filter)
            n_valid = sum(1 for a in record.ALT if a in _VALID_BASES)
            seen_pos[pos0] = seen_pos.get(pos0, 0) + n_valid
            if seen_pos[pos0] > 1:
                multiallelic_pos.add(pos0)

            for alt_idx, alt_allele in enumerate(record.ALT):
                # Per-allele SNP check: single valid nucleotide
                if alt_allele not in _VALID_BASES:
                    self._n_star_alleles += 1
                    continue
                depth, freq, ref_count, alt_count = (
                    self._extract_freq_depth(record, alt_idx))
                # min_depth and min_freq filters are deferred to the
                # grouping phase so that multiallelic merge uses complete
                # allele depths for the denominator.
                raw.append((
                    pos0, record.REF, alt_allele,
                    freq, depth, ref_count, alt_count,
                ))

        if not raw:
            return []

        # Group by position to merge multiallelic decomposed records
        by_pos: dict[int, list[tuple]] = {}
        for rec in raw:
            by_pos.setdefault(rec[0], []).append(rec)

        variants: list[Variant] = []
        for pos in sorted(by_pos):
            group = by_pos[pos]
            is_multiallelic = pos in multiallelic_pos
            if not is_multiallelic:
                p, ref, alt, freq, depth, _, _ = group[0]
                if depth < self._min_depth or freq < self._min_freq:
                    continue
                variants.append(Variant(
                    pos=p, ref=ref, alt=alt, freq=freq, depth=depth))
            elif not self._keep_multiallelic:
                # Skip multiallelic sites by default
                continue
            else:
                # Merge: recompute frequencies using ALL allele depths
                # (including alleles that may fall below min_freq or
                # min_depth) so that the denominator is correct.
                ref_count = group[0][5]
                if ref_count > 0 or any(r[6] > 0 for r in group):
                    total_alt = sum(r[6] for r in group)
                    total_depth = ref_count + total_alt
                    if total_depth < self._min_depth:
                        continue
                    for r in group:
                        new_freq = (r[6] / total_depth
                                    if total_depth > 0 else 0.0)
                        if new_freq >= self._min_freq:
                            variants.append(Variant(
                                pos=r[0], ref=r[1], alt=r[2],
                                freq=new_freq, depth=total_depth))
                else:
                    # No raw allele depths available — keep original frequencies
                    for r in group:
                        if r[4] < self._min_depth or r[3] < self._min_freq:
                            continue
                        variants.append(Variant(
                            pos=r[0], ref=r[1], alt=r[2],
                            freq=r[3], depth=r[4]))

        return variants

    def _extract_freq_depth(
        self, record, alt_index: int = 0,
    ) -> tuple[int, float, int, int]:
        """Extract depth and alt freq for a specific ALT allele.

        Args:
            record: cyvcf2 Variant record.
            alt_index: 0-based index into record.ALT.

        Returns:
            (total_depth, alt_frequency, ref_allele_count, alt_allele_count)
        """
        try:
            ad = record.format("AD")
            if ad is not None:
                ref_count = int(ad[0][0])
                alt_count = int(ad[0][alt_index + 1])
                if ref_count < 0 or alt_count < 0:
                    raise ValueError("missing AD value")
                # Total depth across all alleles at this site
                total = sum(int(x) for x in ad[0] if int(x) >= 0)
                if total > 0:
                    return total, alt_count / total, ref_count, alt_count
        except (KeyError, IndexError, TypeError, ValueError):
            pass

        try:
            dp = record.format("DP")
            depth = int(dp[0][0]) if dp is not None else 0
        except (KeyError, IndexError, TypeError):
            depth = 0
        if depth == 0:
            try:
                depth = record.INFO.get("DP", 0)
            except (KeyError, TypeError):
                depth = 0

        try:
            af = record.INFO.get("AF")
            if isinstance(af, (list, tuple)):
                freq = (float(af[alt_index])
                        if alt_index < len(af) else 0.0)
            else:
                freq = (float(af)
                        if af is not None and alt_index == 0 else 0.0)
        except (KeyError, TypeError, ValueError):
            freq = 0.0

        return int(depth), freq, 0, 0


class IndividualVariantReader(_BaseVariantReader):
    """Variant reader for individual-sequencing data (GT-based frequencies).

    Derives pooled allele frequencies from genotype fields across selected
    diploid samples.  Exposes the same ``fetch()`` interface as
    ``VariantReader`` so downstream code is unchanged.
    """

    def __init__(self, vcf_path: str, samples: list[str] | None = None,
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
        self._n_samples = len(self._vcf.samples)
        self._init_contig_tracking()

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
        if not self._check_contig(chrom):
            return []

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
                    self._n_star_alleles += 1
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
