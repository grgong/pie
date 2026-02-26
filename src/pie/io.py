"""TSV output writers for gene, window, and summary piN/piS results."""

from __future__ import annotations

import numpy as np
import pandas as pd

from pie.diversity import GeneResult


def write_gene_results(results: list[GeneResult], path: str) -> None:
    """Write per-gene TSV with piN, piS, and piN/piS columns."""
    rows = []
    for r in results:
        piN = r.piN
        piS = r.piS
        piN_piS = r.piN_piS  # None when piS == 0
        rows.append({
            "chrom": r.chrom,
            "gene_id": r.gene_id,
            "transcript_id": r.transcript_id,
            "start": r.start,
            "end": r.end,
            "strand": r.strand,
            "n_codons": r.n_codons,
            "n_poly_codons": r.n_poly_codons,
            "N_sites": r.N_sites,
            "S_sites": r.S_sites,
            "N_diffs": r.N_diffs,
            "S_diffs": r.S_diffs,
            "piN": piN,
            "piS": piS,
            "piN_piS": piN_piS if piN_piS is not None else "NA",
            "mean_depth": r.mean_depth,
            "n_variants": r.n_variants,
        })
    df = pd.DataFrame(rows)
    if not df.empty:
        # Convert "NA" strings to actual NaN for proper TSV handling
        df["piN_piS"] = pd.to_numeric(df["piN_piS"], errors="coerce")
    df.to_csv(path, sep="\t", index=False)


def write_window_results(
    results: list[GeneResult],
    path: str,
    window_size: int,
    window_step: int,
) -> None:
    """Slide bp-based windows across each gene's CDS and aggregate codon results."""
    rows = []
    for r in results:
        if not r.codon_results:
            continue
        # Determine CDS span from codon genomic positions
        positions = [cr.pos1 for cr in r.codon_results]
        cds_min = min(positions)
        cds_max = max(positions)

        win_start = cds_min
        while win_start <= cds_max:
            win_end = win_start + window_size
            # Collect codons whose pos1 falls in [win_start, win_end)
            window_codons = [
                cr for cr in r.codon_results
                if win_start <= cr.pos1 < win_end
            ]
            n_codons = len(window_codons)
            N_sites = sum(cr.N_sites for cr in window_codons)
            S_sites = sum(cr.S_sites for cr in window_codons)
            N_diffs = sum(cr.N_diffs for cr in window_codons)
            S_diffs = sum(cr.S_diffs for cr in window_codons)
            piN = N_diffs / N_sites if N_sites > 0 else 0.0
            piS = S_diffs / S_sites if S_sites > 0 else 0.0
            piN_piS = piN / piS if piS > 0 else None

            rows.append({
                "chrom": r.chrom,
                "win_start": win_start,
                "win_end": win_end,
                "gene_id": r.gene_id,
                "n_codons": n_codons,
                "N_sites": N_sites,
                "S_sites": S_sites,
                "N_diffs": N_diffs,
                "S_diffs": S_diffs,
                "piN": piN,
                "piS": piS,
                "piN_piS": piN_piS if piN_piS is not None else "NA",
            })
            win_start += window_step

    df = pd.DataFrame(rows)
    if not df.empty:
        df["piN_piS"] = pd.to_numeric(df["piN_piS"], errors="coerce")
    df.to_csv(path, sep="\t", index=False)


def write_summary(results: list[GeneResult], path: str) -> None:
    """Write single-row genome-wide summary TSV."""
    total_genes = len(results)
    total_codons = sum(r.n_codons for r in results)
    total_variants = sum(r.n_variants for r in results)

    total_N_sites = sum(r.N_sites for r in results)
    total_S_sites = sum(r.S_sites for r in results)
    total_N_diffs = sum(r.N_diffs for r in results)
    total_S_diffs = sum(r.S_diffs for r in results)

    genome_piN = total_N_diffs / total_N_sites if total_N_sites > 0 else 0.0
    genome_piS = total_S_diffs / total_S_sites if total_S_sites > 0 else 0.0
    genome_piN_piS = genome_piN / genome_piS if genome_piS > 0 else None

    gene_piNs = [r.piN for r in results]
    gene_piSs = [r.piS for r in results]

    row = {
        "total_genes": total_genes,
        "total_codons": total_codons,
        "total_variants": total_variants,
        "genome_piN": genome_piN,
        "genome_piS": genome_piS,
        "genome_piN_piS": genome_piN_piS if genome_piN_piS is not None else "NA",
        "mean_gene_piN": float(np.mean(gene_piNs)) if gene_piNs else 0.0,
        "mean_gene_piS": float(np.mean(gene_piSs)) if gene_piSs else 0.0,
        "median_gene_piN": float(np.median(gene_piNs)) if gene_piNs else 0.0,
        "median_gene_piS": float(np.median(gene_piSs)) if gene_piSs else 0.0,
    }
    df = pd.DataFrame([row])
    df["genome_piN_piS"] = pd.to_numeric(df["genome_piN_piS"], errors="coerce")
    df.to_csv(path, sep="\t", index=False)
