import pytest
from pie.reference import ReferenceGenome


class TestReferenceGenome:
    def test_context_manager(self, ref_fasta):
        with ReferenceGenome(ref_fasta) as ref:
            assert ref is not None

    def test_fetch_sequence(self, ref_fasta):
        with ReferenceGenome(ref_fasta) as ref:
            seq = ref.fetch("chr1", 0, 3)
            assert seq == "ATG"  # Gene1 starts with ATG

    def test_extract_codons_plus_strand(self, ref_fasta):
        with ReferenceGenome(ref_fasta) as ref:
            exons = [("chr1", 0, 90)]
            codons = ref.extract_codons(exons, strand="+")
            assert len(codons) == 30
            assert codons[0] == "ATG"
            assert codons[1] == "GCT"
            assert codons[2] == "GAT"
            assert codons[-1] == "TAA"

    def test_extract_codons_minus_strand(self, ref_fasta):
        with ReferenceGenome(ref_fasta) as ref:
            exons = [("chr1", 230, 311)]
            codons = ref.extract_codons(exons, strand="-")
            assert len(codons) == 27
            assert codons[0] == "ATG"  # First codon should be start

    def test_extract_codons_multi_exon(self, ref_fasta):
        with ReferenceGenome(ref_fasta) as ref:
            exons = [("chr1", 100, 159), ("chr1", 180, 220)]
            codons = ref.extract_codons(exons, strand="+")
            assert len(codons) == 33  # 99bp / 3
            assert codons[0] == "ATG"
            # Exon-spanning codon (codon 20): last 2bp of exon1 + first 1bp of exon2
            assert codons[19] == "ACT"  # Thr

    def test_codon_genomic_positions_plus(self, ref_fasta):
        with ReferenceGenome(ref_fasta) as ref:
            exons = [("chr1", 0, 90)]
            positions = ref.codon_genomic_positions(exons, strand="+")
            assert len(positions) == 30
            assert positions[0] == ("chr1", 0, 1, 2)
            assert positions[1] == ("chr1", 3, 4, 5)
            assert positions[29] == ("chr1", 87, 88, 89)

    def test_codon_genomic_positions_minus(self, ref_fasta):
        with ReferenceGenome(ref_fasta) as ref:
            exons = [("chr1", 230, 311)]
            positions = ref.codon_genomic_positions(exons, strand="-")
            assert len(positions) == 27
            # First codon on - strand: rightmost 3 bases reversed
            assert positions[0] == ("chr1", 310, 309, 308)

    def test_codon_genomic_positions_multi_exon(self, ref_fasta):
        with ReferenceGenome(ref_fasta) as ref:
            exons = [("chr1", 100, 159), ("chr1", 180, 220)]
            positions = ref.codon_genomic_positions(exons, strand="+")
            assert len(positions) == 33
            # Exon-spanning codon 20 (0-indexed 19):
            # last 2bp of exon1 (157, 158) + first bp of exon2 (180)
            assert positions[19] == ("chr1", 157, 158, 180)
