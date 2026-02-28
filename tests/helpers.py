"""Shared test utilities for creating bgzipped VCFs and indexed FASTAs."""

import subprocess


def bgzip_and_index(vcf_path):
    """Bgzip and tabix-index a plain VCF file. Returns .vcf.gz path."""
    gz_path = str(vcf_path) + ".gz"
    with open(gz_path, "wb") as out:
        subprocess.run(["bgzip", "-c", str(vcf_path)], stdout=out, check=True)
    subprocess.run(["tabix", "-p", "vcf", gz_path], check=True)
    return gz_path


def write_fasta(path, sequences: dict[str, str]):
    """Write a multi-sequence FASTA and samtools-index it."""
    with open(path, "w") as fh:
        for name, seq in sequences.items():
            fh.write(f">{name}\n{seq}\n")
    subprocess.run(["samtools", "faidx", str(path)], check=True)
    return str(path)
