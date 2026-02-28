"""Reference genome FASTA access and codon extraction."""

import pysam

_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


class ReferenceGenome:
    def __init__(self, fasta_path: str):
        self._fasta = pysam.FastaFile(fasta_path)

    def close(self):
        self._fasta.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    def fetch(self, chrom: str, start: int, end: int) -> str:
        """Fetch sequence (0-based, half-open)."""
        return self._fasta.fetch(chrom, start, end).upper()

    def extract_codons(
        self, exons: list[tuple[str, int, int]], strand: str
    ) -> list[str]:
        """Extract codons from CDS exons.

        Args:
            exons: [(chrom, start, end), ...] in genomic order, 0-based half-open.
            strand: "+" or "-"
        Returns:
            List of 3-letter codon strings in reading frame order.
        """
        cds_seq = "".join(self.fetch(chrom, start, end) for chrom, start, end in exons)
        if strand == "-":
            cds_seq = cds_seq[::-1].translate(_COMPLEMENT)
        n_complete = (len(cds_seq) // 3) * 3
        return [cds_seq[i : i + 3] for i in range(0, n_complete, 3)]

    def codon_genomic_positions(
        self, exons: list[tuple[str, int, int]], strand: str
    ) -> list[tuple[str, int, int, int]]:
        """Map codon indices to genomic positions.

        Returns: [(chrom, pos1, pos2, pos3), ...] where pos are 0-based genomic.
        """
        chrom_list: list[str] = []
        pos_list: list[int] = []
        for chrom, start, end in exons:
            n = end - start
            chrom_list.extend([chrom] * n)
            pos_list.extend(range(start, end))
        if strand == "-":
            chrom_list.reverse()
            pos_list.reverse()
        n_complete = (len(pos_list) // 3) * 3
        return [
            (chrom_list[i], pos_list[i], pos_list[i + 1], pos_list[i + 2])
            for i in range(0, n_complete, 3)
        ]
