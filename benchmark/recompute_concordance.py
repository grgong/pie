#!/usr/bin/env python
"""Recompute concordance metrics between pie and SNPGenie outputs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd
from scipy.stats import pearsonr, spearmanr


def _read_snpgenie_products(root: Path, rep: int, chroms: list[str]) -> pd.DataFrame:
    frames = []
    for chrom in chroms:
        path = root / f"{chrom}_rep{rep}" / "SNPGenie_Results" / "product_results.txt"
        df = pd.read_csv(path, sep="\t")
        df["chrom"] = chrom
        frames.append(df)
    out = pd.concat(frames, ignore_index=True)
    for col in ["N_diffs", "S_diffs", "N_sites", "S_sites", "piN", "piS"]:
        out[col] = pd.to_numeric(out[col], errors="coerce")
    return out


def _read_snpgenie_codon_rows(root: Path, rep: int, chrom: str, product: str) -> pd.DataFrame:
    path = root / f"{chrom}_rep{rep}" / "SNPGenie_Results" / "codon_results.txt"
    df = pd.read_csv(path, sep="\t")
    return df[df["product"] == product].copy()


def _read_pie(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    df["product"] = df["gene_id"].str.replace(r"^gene-", "", regex=True)
    return df


def _weighted_metric(df: pd.DataFrame, diffs_col: str, sites_col: str) -> float:
    return float(df[diffs_col].sum() / df[sites_col].sum())


def _rel_diff_percent(a: float, b: float) -> float:
    return abs(a - b) / b * 100.0


def _corr(x: pd.Series, y: pd.Series) -> tuple[float, float]:
    return float(pearsonr(x, y).statistic), float(spearmanr(x, y).statistic)


def _percentage(mask: pd.Series) -> float:
    return float(mask.mean() * 100.0)


def _outlier_vcf_records(vcf_path: Path, chrom: str, start: int, end: int) -> list[dict[str, object]]:
    records = []
    with open(vcf_path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            rec_chrom, pos = fields[0], int(fields[1])
            if rec_chrom != chrom or pos < start or pos > end:
                continue
            fmt_keys = fields[8].split(":")
            sample_vals = fields[9].split(":")
            fmt = dict(zip(fmt_keys, sample_vals))
            records.append(
                {
                    "pos": pos,
                    "ref": fields[3],
                    "alt": fields[4],
                    "ro": int(fmt.get("RO", "0") or 0),
                    "ao": int((fmt.get("AO", "0") or "0").split(",")[0]),
                    "dp": int(fmt.get("DP", "0") or 0),
                }
            )
    return records


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pie-gene-results", required=True, type=Path)
    parser.add_argument("--snpgenie-root", required=True, type=Path)
    parser.add_argument("--rep", default=1, type=int)
    parser.add_argument("--pie-gene-results-no-keep", type=Path, default=None)
    parser.add_argument("--vcf", type=Path, default=None)
    parser.add_argument("--output-json", type=Path, default=None)
    parser.add_argument("--chroms", nargs="+", default=["X", "1", "2", "3"])
    args = parser.parse_args()

    pie = _read_pie(args.pie_gene_results)
    snp = _read_snpgenie_products(args.snpgenie_root, args.rep, args.chroms)

    merged = pie.merge(
        snp[["product", "chrom", "N_diffs", "S_diffs", "N_sites", "S_sites", "piN", "piS"]],
        on="product",
        how="inner",
        suffixes=("_pie", "_snp"),
    )

    genome_pin_pie = _weighted_metric(merged, "N_diffs_pie", "N_sites_pie")
    genome_pis_pie = _weighted_metric(merged, "S_diffs_pie", "S_sites_pie")
    genome_pin_snp = _weighted_metric(merged, "N_diffs_snp", "N_sites_snp")
    genome_pis_snp = _weighted_metric(merged, "S_diffs_snp", "S_sites_snp")

    piN_rel = (
        (merged.loc[(merged["piN_pie"] > 0) & (merged["piN_snp"] > 0), "piN_pie"]
         - merged.loc[(merged["piN_pie"] > 0) & (merged["piN_snp"] > 0), "piN_snp"])
        .abs()
        / merged.loc[(merged["piN_pie"] > 0) & (merged["piN_snp"] > 0), "piN_snp"]
    )
    piS_rel = (
        (merged.loc[(merged["piS_pie"] > 0) & (merged["piS_snp"] > 0), "piS_pie"]
         - merged.loc[(merged["piS_pie"] > 0) & (merged["piS_snp"] > 0), "piS_snp"])
        .abs()
        / merged.loc[(merged["piS_pie"] > 0) & (merged["piS_snp"] > 0), "piS_snp"]
    )

    outlier = merged.loc[((merged["piS_pie"] > 0) & (merged["piS_snp"] > 0))].copy()
    known_outlier = outlier[outlier["product"] == "Apisum_003668"]
    if not known_outlier.empty:
        outlier_row = known_outlier.iloc[0]
    else:
        outlier["piS_rel_diff"] = (outlier["piS_pie"] - outlier["piS_snp"]).abs() / outlier["piS_snp"]
        outlier_row = outlier.sort_values("piS_rel_diff", ascending=False).iloc[0]
    outlier_product = str(outlier_row["product"])
    outlier_chrom = str(outlier_row["chrom_pie"])

    piS_ex_outlier = merged[merged["product"] != outlier_product]
    outlier_codons = _read_snpgenie_codon_rows(args.snpgenie_root, args.rep, outlier_chrom, outlier_product)
    outlier_codon = outlier_codons.sort_values("S_diffs", ascending=False).iloc[0]

    outlier_vcf = []
    if args.vcf is not None:
        outlier_vcf = _outlier_vcf_records(args.vcf, outlier_chrom, int(outlier_codon["site"]), int(outlier_codon["site"]) + 2)

    metrics = {
        "matched_genes": int(len(merged)),
        "genome": {
            "pie": {
                "piN": genome_pin_pie,
                "piS": genome_pis_pie,
                "piN_piS": genome_pin_pie / genome_pis_pie,
            },
            "snpgenie": {
                "piN": genome_pin_snp,
                "piS": genome_pis_snp,
                "piN_piS": genome_pin_snp / genome_pis_snp,
            },
            "relative_diff_percent": {
                "piN": _rel_diff_percent(genome_pin_pie, genome_pin_snp),
                "piS": _rel_diff_percent(genome_pis_pie, genome_pis_snp),
                "piN_piS": _rel_diff_percent(genome_pin_pie / genome_pis_pie, genome_pin_snp / genome_pis_snp),
            },
        },
        "per_gene": {
            "site_counts": {
                "N_sites": {
                    "pearson": _corr(merged["N_sites_pie"], merged["N_sites_snp"])[0],
                    "spearman": _corr(merged["N_sites_pie"], merged["N_sites_snp"])[1],
                },
                "S_sites": {
                    "pearson": _corr(merged["S_sites_pie"], merged["S_sites_snp"])[0],
                    "spearman": _corr(merged["S_sites_pie"], merged["S_sites_snp"])[1],
                },
            },
            "piN": {
                "pearson": _corr(merged["piN_pie"], merged["piN_snp"])[0],
                "spearman": _corr(merged["piN_pie"], merged["piN_snp"])[1],
            },
            "piS": {
                "pearson_excluding_outlier": _corr(piS_ex_outlier["piS_pie"], piS_ex_outlier["piS_snp"])[0],
                "spearman_excluding_outlier": _corr(piS_ex_outlier["piS_pie"], piS_ex_outlier["piS_snp"])[1],
            },
            "N_diffs_S_diffs": {
                "pearson": _corr(merged["N_diffs_pie"], merged["N_diffs_snp"])[0],
                "spearman": _corr(merged["N_diffs_pie"], merged["N_diffs_snp"])[1],
            },
            "median_relative_diff_percent": {
                "piN": float(piN_rel.median() * 100.0),
                "piS": float(piS_rel.median() * 100.0),
            },
            "within_5_percent": {
                "piN": _percentage(piN_rel <= 0.05),
                "piS": _percentage(piS_rel <= 0.05),
            },
        },
        "outlier_gene": {
            "product": outlier_product,
            "chrom": outlier_chrom,
            "pie_piS": float(outlier_row["piS_pie"]),
            "snpgenie_piS": float(outlier_row["piS_snp"]),
            "piS_fold_difference": float(outlier_row["piS_snp"] / outlier_row["piS_pie"]),
            "snpgenie_max_codon": {
                "site": int(outlier_codon["site"]),
                "codon": str(outlier_codon["codon"]),
                "S_diffs": float(outlier_codon["S_diffs"]),
                "N_sites": float(outlier_codon["N_sites"]),
                "S_sites": float(outlier_codon["S_sites"]),
            },
            "vcf_records_at_outlier_codon": outlier_vcf,
            "vcf_records_with_zero_ref": int(sum(rec["ro"] == 0 for rec in outlier_vcf)),
        },
    }

    if args.pie_gene_results_no_keep is not None:
        pie_no_keep = _read_pie(args.pie_gene_results_no_keep)
        merged_no_keep = pie_no_keep.merge(
            snp[["product", "N_diffs", "S_diffs", "N_sites", "S_sites"]],
            on="product",
            how="inner",
            suffixes=("_pie", "_snp"),
        )
        genome_pin_no = _weighted_metric(merged_no_keep, "N_diffs_pie", "N_sites_pie")
        metrics["multiallelic_effect"] = {
            "genes_with_filtered_multiallelic": int((pie_no_keep["n_filtered_multiallelic"] > 0).sum()),
            "filtered_multiallelic_records": int(pie_no_keep["n_filtered_multiallelic"].sum()),
            "matched_gene_piN_without_keep": genome_pin_no,
            "matched_gene_piN_with_keep": genome_pin_pie,
            "relative_diff_vs_snpgenie_without_keep_percent": _rel_diff_percent(genome_pin_no, genome_pin_snp),
            "relative_diff_vs_snpgenie_with_keep_percent": _rel_diff_percent(genome_pin_pie, genome_pin_snp),
        }

    if args.output_json is not None:
        args.output_json.write_text(json.dumps(metrics, indent=2) + "\n")

    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
