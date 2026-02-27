#!/usr/bin/env python3
"""Convert GFF3 CDS features to SNPGenie-compatible GTF format.

SNPGenie requires GTF with gene_id and transcript_id attributes.
We extract locus_tag from GFF3 as gene_id and orig_transcript_id as transcript_id.
Only one isoform per gene is kept (the one with the most CDS bp).
"""
import sys
import re
from collections import defaultdict


def parse_gff3_attributes(attr_str):
    attrs = {}
    for item in attr_str.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            attrs[k.strip()] = v.strip()
    return attrs


def main(gff3_path, gtf_path):
    # Collect CDS lines grouped by transcript
    transcripts = defaultdict(list)  # transcript_id -> [(chrom, src, start, end, score, strand, frame, gene_id)]

    with open(gff3_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "CDS":
                continue
            attrs = parse_gff3_attributes(parts[8])
            gene_id = attrs.get("locus_tag", "")
            transcript_id = attrs.get("orig_transcript_id", attrs.get("Parent", ""))
            # Clean transcript_id: remove "rna-" prefix if present
            transcript_id = re.sub(r"^rna-", "", transcript_id)
            if not gene_id:
                continue
            transcripts[transcript_id].append(
                (parts[0], parts[1], int(parts[3]), int(parts[4]),
                 parts[5], parts[6], parts[7], gene_id)
            )

    # Pick longest isoform per gene
    gene_isoforms = defaultdict(list)  # gene_id -> [(transcript_id, total_bp, cds_list)]
    for tid, cds_list in transcripts.items():
        gene_id = cds_list[0][7]
        total_bp = sum(end - start + 1 for _, _, start, end, _, _, _, _ in cds_list)
        gene_isoforms[gene_id].append((tid, total_bp, cds_list))

    with open(gtf_path, "w") as out:
        for gene_id, isoforms in sorted(gene_isoforms.items()):
            # Pick longest
            best = max(isoforms, key=lambda x: x[1])
            tid, _, cds_list = best
            for chrom, src, start, end, score, strand, frame, _ in sorted(cds_list, key=lambda x: x[2]):
                attr = f'gene_id "{gene_id}"; transcript_id "{tid}";'
                out.write(f"{chrom}\t{src}\tCDS\t{start}\t{end}\t{score}\t{strand}\t{frame}\t{attr}\n")

    n_genes = len(gene_isoforms)
    n_cds = sum(1 for isos in gene_isoforms.values() for _ in max(isos, key=lambda x: x[1])[2])
    print(f"Converted {n_cds} CDS features from {n_genes} genes to GTF")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
