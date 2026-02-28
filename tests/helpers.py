"""Shared test utilities for creating bgzipped VCFs and indexed FASTAs."""

import subprocess

import pysam


def bgzip_and_index(vcf_path):
    """Bgzip and tabix-index a plain VCF file. Returns .vcf.gz path."""
    gz_path = str(vcf_path) + ".gz"
    pysam.tabix_compress(str(vcf_path), gz_path, force=True)
    pysam.tabix_index(gz_path, preset="vcf", force=True)
    return gz_path


def write_fasta(path, sequences: dict[str, str]):
    """Write a multi-sequence FASTA and samtools-index it."""
    with open(path, "w") as fh:
        for name, seq in sequences.items():
            fh.write(f">{name}\n{seq}\n")
    subprocess.run(["samtools", "faidx", str(path)], check=True)
    return str(path)
