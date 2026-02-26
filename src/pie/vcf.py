"""VCF parsing, filtering, and auto-indexing via cyvcf2."""
import os
import subprocess
import logging
from dataclasses import dataclass
from cyvcf2 import VCF

log = logging.getLogger(__name__)

_VALID_BASES = frozenset("ACGT")


@dataclass(slots=True)
class Variant:
    pos: int      # 0-based genomic position
    ref: str      # reference allele
    alt: str      # alternate allele
    freq: float   # alt allele frequency (0-1)
    depth: int    # total depth


def ensure_indexed(vcf_path: str) -> str:
    """Ensure VCF is bgzipped and tabix-indexed. Returns path to .vcf.gz."""
    if vcf_path.endswith(".vcf"):
        gz_path = vcf_path + ".gz"
        log.info("Bgzipping %s -> %s", vcf_path, gz_path)
        with open(gz_path, "wb") as out:
            subprocess.run(["bgzip", "-c", vcf_path], stdout=out, check=True)
        vcf_path = gz_path

    tbi_path = vcf_path + ".tbi"
    if not os.path.exists(tbi_path):
        log.info("Creating tabix index for %s", vcf_path)
        subprocess.run(["tabix", "-p", "vcf", vcf_path], check=True)

    return vcf_path


class VariantReader:
    def __init__(self, vcf_path: str, min_freq: float = 0.01,
                 min_depth: int = 10, min_qual: float = 20.0,
                 pass_only: bool = False, keep_multiallelic: bool = False):
        self._vcf_path = vcf_path
        self._min_freq = min_freq
        self._min_depth = min_depth
        self._min_qual = min_qual
        self._pass_only = pass_only
        self._keep_multiallelic = keep_multiallelic
        self._vcf = VCF(vcf_path)

    def close(self):
        self._vcf.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    def fetch(self, chrom: str, start: int, end: int) -> list[Variant]:
        """Fetch filtered variants in region (0-based, half-open).

        Handles both single-line multiallelic records (``ALT=C,G``) and
        decomposed records (from ``bcftools norm -m-``).  By default,
        positions with multiple ALT alleles are skipped; set
        ``keep_multiallelic=True`` to merge them instead (frequencies
        recomputed from raw allele depths).
        """
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
                    continue
                depth, freq, ref_count, alt_count = (
                    self._extract_freq_depth(record, alt_idx))
                if depth < self._min_depth:
                    continue
                if freq < self._min_freq:
                    continue
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
                variants.append(Variant(
                    pos=p, ref=ref, alt=alt, freq=freq, depth=depth))
            elif not self._keep_multiallelic:
                # Skip multiallelic sites by default
                continue
            else:
                # Merge: recompute frequencies using allele depths
                ref_count = group[0][5]
                if ref_count > 0 or any(r[6] > 0 for r in group):
                    total_alt = sum(r[6] for r in group)
                    total_depth = ref_count + total_alt
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
