"""VCF parsing, filtering, and auto-indexing via cyvcf2."""
import os
import subprocess
import logging
from dataclasses import dataclass
from cyvcf2 import VCF

log = logging.getLogger(__name__)


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
                 pass_only: bool = False):
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
        """Fetch filtered variants in region (0-based, half-open)."""
        variants = []
        region = f"{chrom}:{start + 1}-{end}"
        for record in self._vcf(region):
            if not record.is_snp:
                continue
            if len(record.ALT) != 1:
                continue
            if self._pass_only and record.FILTER is not None:
                continue
            if record.QUAL is not None and record.QUAL < self._min_qual:
                continue
            depth, freq = self._extract_freq_depth(record)
            if depth < self._min_depth:
                continue
            if freq < self._min_freq:
                continue
            variants.append(Variant(
                pos=record.POS - 1, ref=record.REF,
                alt=record.ALT[0], freq=freq, depth=depth,
            ))
        return variants

    def _extract_freq_depth(self, record) -> tuple[int, float]:
        """Extract depth and alt freq. Tries AD first, then INFO AF/DP."""
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
                freq = float(af[0])
            else:
                freq = float(af) if af is not None else 0.0
        except (KeyError, TypeError, ValueError):
            freq = 0.0

        return int(depth), freq
