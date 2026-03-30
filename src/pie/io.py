"""TSV output writers for gene, window, and summary piN/piS results."""

from __future__ import annotations

import dataclasses
from bisect import bisect_left

import numpy as np
import pandas as pd

from pie.diversity import GeneResult, VariantRecord

_VARIANT_COLUMNS = [f.name for f in dataclasses.fields(VariantRecord)]


def write_variant_results(records: list[VariantRecord], path: str) -> None:
    """Write per-variant TSV with codon annotation and allele counts."""
    if not records:
        pd.DataFrame(columns=_VARIANT_COLUMNS).to_csv(path, sep="\t", index=False)
        return
    rows = [[getattr(r, c) for c in _VARIANT_COLUMNS] for r in records]
    pd.DataFrame(rows, columns=_VARIANT_COLUMNS).to_csv(path, sep="\t", index=False)


def write_gene_results(results: list[GeneResult], path: str) -> None:
    """Write per-gene TSV with piN, piS, and piN/piS columns."""
    ind_mode = any(r.n_samples is not None for r in results)
    rows = []
    for r in results:
        piN = r.piN
        piS = r.piS
        piN_piS = r.piN_piS  # None when piS == 0
        fs = r.filter_stats
        row = {
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
            "piN_piS": piN_piS,
            "mean_variant_depth": r.mean_variant_depth,
            "n_variants": r.n_variants,
            # QC columns
            "n_ambiguous_codons": r.n_ambiguous_codons,
            "n_internal_stop_codons": r.n_internal_stop_codons,
            "n_vcf_records": fs.n_total,
            "n_filtered_qual": fs.n_filtered_qual,
            "n_filtered_depth": fs.n_filtered_depth,
            "n_filtered_freq": fs.n_filtered_freq,
            "n_filtered_multiallelic": fs.n_filtered_multiallelic,
            "n_filtered_not_snp": fs.n_filtered_not_snp,
        }
        if ind_mode:
            row["n_samples"] = r.n_samples
            row["mean_call_rate"] = r.mean_call_rate
            row["n_filtered_call_rate"] = fs.n_filtered_call_rate
        rows.append(row)
    df = pd.DataFrame(rows)
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

        # Build sorted position array and prefix sums for O(log n) window queries
        codons = sorted(r.codon_results, key=lambda cr: cr.pos1)
        n = len(codons)
        pos_arr = [cr.pos1 for cr in codons]
        cds_min = pos_arr[0]
        cds_max = pos_arr[-1]

        # Prefix sums via numpy (length n+1, prefix[0] = 0)
        vals = np.array(
            [(cr.N_sites, cr.S_sites, cr.N_diffs, cr.S_diffs) for cr in codons],
            dtype=np.float64,
        )
        prefix = np.zeros((n + 1, 4), dtype=np.float64)
        np.cumsum(vals, axis=0, out=prefix[1:])

        chrom = r.chrom
        gene_id = r.gene_id
        win_start = cds_min
        while win_start <= cds_max:
            win_end = win_start + window_size
            # Bisect for codons with pos1 in [win_start, win_end)
            lo = bisect_left(pos_arr, win_start)
            hi = bisect_left(pos_arr, win_end)
            n_codons = hi - lo
            N_sites = prefix[hi, 0] - prefix[lo, 0]
            S_sites = prefix[hi, 1] - prefix[lo, 1]
            N_diffs = prefix[hi, 2] - prefix[lo, 2]
            S_diffs = prefix[hi, 3] - prefix[lo, 3]
            piN = N_diffs / N_sites if N_sites > 0 else 0.0
            piS = S_diffs / S_sites if S_sites > 0 else 0.0
            piN_piS = piN / piS if piS > 0 else None

            rows.append({
                "chrom": chrom,
                "win_start": win_start,
                "win_end": win_end,
                "gene_id": gene_id,
                "n_codons": n_codons,
                "N_sites": N_sites,
                "S_sites": S_sites,
                "N_diffs": N_diffs,
                "S_diffs": S_diffs,
                "piN": piN,
                "piS": piS,
                "piN_piS": piN_piS,
            })
            win_start += window_step

    df = pd.DataFrame(rows)
    df.to_csv(path, sep="\t", index=False)


def write_summary(results: list[GeneResult], path: str) -> None:
    """Write single-row genome-wide summary TSV."""
    total_genes = len(results)
    total_codons = sum(r.n_codons for r in results)
    cds_snp_variants = sum(r.n_variants for r in results)

    total_N_sites = sum(r.N_sites for r in results)
    total_S_sites = sum(r.S_sites for r in results)
    total_N_diffs = sum(r.N_diffs for r in results)
    total_S_diffs = sum(r.S_diffs for r in results)

    genome_piN = total_N_diffs / total_N_sites if total_N_sites > 0 else 0.0
    genome_piS = total_S_diffs / total_S_sites if total_S_sites > 0 else 0.0
    genome_piN_piS = genome_piN / genome_piS if genome_piS > 0 else None

    gene_piNs = [r.piN for r in results]
    gene_piSs = [r.piS for r in results]

    ind_mode = any(r.n_samples is not None for r in results)

    total_vcf_records = sum(r.filter_stats.n_total for r in results)
    total_filtered_qual = sum(r.filter_stats.n_filtered_qual for r in results)
    total_filtered_depth = sum(r.filter_stats.n_filtered_depth for r in results)
    total_filtered_freq = sum(r.filter_stats.n_filtered_freq for r in results)
    total_filtered_multiallelic = sum(
        r.filter_stats.n_filtered_multiallelic for r in results)
    total_filtered_not_snp = sum(
        r.filter_stats.n_filtered_not_snp for r in results)
    total_ambiguous_codons = sum(r.n_ambiguous_codons for r in results)
    total_internal_stop_codons = sum(r.n_internal_stop_codons for r in results)

    row = {
        "total_genes": total_genes,
        "total_codons": total_codons,
        "cds_snp_variants": cds_snp_variants,
        "genome_piN": genome_piN,
        "genome_piS": genome_piS,
        "genome_piN_piS": genome_piN_piS,
        "mean_gene_piN": float(np.mean(gene_piNs)) if gene_piNs else 0.0,
        "mean_gene_piS": float(np.mean(gene_piSs)) if gene_piSs else 0.0,
        "median_gene_piN": float(np.median(gene_piNs)) if gene_piNs else 0.0,
        "median_gene_piS": float(np.median(gene_piSs)) if gene_piSs else 0.0,
        "stop_renorm_genes": sum(1 for r in results if r.n_stop_codons > 0),
        "stop_renorm_codons": sum(r.n_stop_codons for r in results),
        # QC totals
        "total_vcf_records": total_vcf_records,
        "total_ambiguous_codons": total_ambiguous_codons,
        "total_internal_stop_codons": total_internal_stop_codons,
        "total_filtered_qual": total_filtered_qual,
        "total_filtered_depth": total_filtered_depth,
        "total_filtered_freq": total_filtered_freq,
        "total_filtered_multiallelic": total_filtered_multiallelic,
        "total_filtered_not_snp": total_filtered_not_snp,
    }
    if ind_mode:
        row["n_samples_selected"] = results[0].n_samples if results else 0
        all_call_rates = []
        for r in results:
            if r.call_rates:
                all_call_rates.extend(r.call_rates)
        row["mean_call_rate"] = (
            sum(all_call_rates) / len(all_call_rates) if all_call_rates else 0.0
        )
        row["total_filtered_call_rate"] = sum(
            r.filter_stats.n_filtered_call_rate for r in results)
    df = pd.DataFrame([row])
    df.to_csv(path, sep="\t", index=False)
