#!/usr/bin/env python3
"""Convert GFF3 CDS features to SNPGenie-compatible GTF format.

Uses pie's annotation parser to select the longest isoform per gene,
ensuring SNPGenie receives the exact same gene models that pie uses.
"""
import sys

from pie.annotation import parse_annotations


def main(gff3_path, gtf_path):
    genes = parse_annotations(gff3_path)

    n_cds = 0
    with open(gtf_path, "w") as out:
        for gene in genes:
            gene_id = gene.gene_id.replace("gene-", "")
            transcript_id = gene.transcript_id.replace("rna-", "")
            for chrom, start, end in gene.cds_exons:
                # GTF is 1-based inclusive
                attr = f'gene_id "{gene_id}"; transcript_id "{transcript_id}";'
                out.write(
                    f"{chrom}\tpie\tCDS\t{start + 1}\t{end}\t.\t{gene.strand}\t.\t{attr}\n"
                )
                n_cds += 1

    print(f"Converted {n_cds} CDS features from {len(genes)} genes to GTF")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <input.gff3> <output.gtf>", file=sys.stderr)
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
