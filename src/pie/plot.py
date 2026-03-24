"""Publication-ready plots for per-gene and sliding-window piN/piS results."""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
from plotnine import (
    aes,
    after_stat,
    element_blank,
    element_text,
    geom_abline,
    geom_boxplot,
    geom_density,
    geom_histogram,
    geom_hline,
    geom_line,
    geom_point,
    geom_text,
    geom_vline,
    ggplot,
    ggsave,
    guides,
    facet_wrap,
    labs,
    scale_color_manual,
    scale_fill_manual,
    scale_x_continuous,
    theme,
    theme_bw,
)

_OKABE_ITO = [
    "#E69F00", "#56B4E9", "#009E73", "#F0E442",
    "#0072B2", "#D55E00", "#CC79A7", "#999999",
]

# NPG (Nature Publishing Group) figure widths in inches.
NPG_SINGLE = 3.5    # 89 mm — single column
NPG_ONE_HALF = 4.7  # 120 mm — 1.5 column
NPG_DOUBLE = 7.2    # 183 mm — double column


def _base_theme(width: float = NPG_DOUBLE, height: float = 4.0, dpi: int = 300) -> theme:
    """Publication-quality base theme sized for NPG journals."""
    return (
        theme_bw()
        + theme(
            figure_size=(width, height),
            dpi=dpi,
            axis_text=element_text(size=7),
            axis_title=element_text(size=8),
            plot_title=element_text(size=9, weight="bold"),
            panel_grid_minor=element_blank(),
        )
    )


def _sort_chroms(chroms) -> list[str]:
    """Sort chromosome names: numeric first (chr1, chr2, ...), then alphabetic."""
    return sorted(
        chroms,
        key=lambda c: (
            not str(c).replace("chr", "").isdigit(),
            int(str(c).replace("chr", "")) if str(c).replace("chr", "").isdigit() else str(c),
        ),
    )


def _require_columns(df: pd.DataFrame, required: list[str], path: str) -> None:
    """Raise ValueError if any required columns are missing."""
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(
            f"Input file missing required columns: {', '.join(missing)} ({path})"
        )


def _load_gene_results(path: str) -> pd.DataFrame:
    """Load gene_results.tsv and prepare for plotting."""
    df = pd.read_csv(path, sep="\t")
    _require_columns(df, ["chrom", "gene_id", "start", "end"], path)
    df["chrom"] = df["chrom"].astype(str)
    return df


def _apply_base_qc(
    df: pd.DataFrame,
    min_codons: int | None = None,
    min_variants: int | None = None,
) -> pd.DataFrame:
    """Apply shared reliability filters (gene length, variant count)."""
    if min_codons is not None and "n_codons" in df.columns:
        df = df[df["n_codons"] >= min_codons]
    if min_variants is not None and "n_variants" in df.columns:
        df = df[df["n_variants"] >= min_variants]
    return df


def _filter_ratio(
    df: pd.DataFrame,
    max_ratio: float = 2.0,
    exclude_zero_ratio: bool = False,
) -> pd.DataFrame:
    """Filter for piN/piS analysis: !is.na & piS>0 & piN_piS>=0 & piN_piS<=max."""
    df = df.dropna(subset=["piN_piS"])
    if "piS" in df.columns:
        df = df[df["piS"] > 0]
    df = df[(df["piN_piS"] >= 0) & (df["piN_piS"] <= max_ratio)]
    if exclude_zero_ratio:
        df = df[df["piN_piS"] > 0]
    return df


def _filter_metric(
    df: pd.DataFrame,
    column: str,
    max_value: float = 2.0,
) -> pd.DataFrame:
    """Filter a single metric: !is.na & value<=max."""
    df = df.dropna(subset=[column])
    df = df[df[column] <= max_value]
    return df


def _add_cumulative_pos(df: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, float]]:
    """Add cumulative genomic position column and return chrom center positions."""
    chroms = _sort_chroms(df["chrom"].unique())
    chrom_offsets = {}
    chrom_centers = {}
    cumulative = 0
    for chrom in chroms:
        chrom_offsets[chrom] = cumulative
        chrom_max = df.loc[df["chrom"] == chrom, "end"].max()
        chrom_centers[chrom] = cumulative + chrom_max / 2
        cumulative += chrom_max * 1.05

    df = df.copy()
    df["cum_pos"] = df.apply(
        lambda row: chrom_offsets[row["chrom"]] + (row["start"] + row["end"]) / 2,
        axis=1,
    )
    chrom_to_idx = {c: i for i, c in enumerate(chroms)}
    df["chrom_color"] = df["chrom"].map(lambda c: _OKABE_ITO[chrom_to_idx[c] % len(_OKABE_ITO)])
    df["chrom"] = pd.Categorical(df["chrom"], categories=chroms, ordered=True)
    return df, chrom_centers


def _save_plot(p, output_path: str, dpi: int) -> None:
    """Save a plotnine plot, suppressing warnings."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ggsave(p, output_path, dpi=dpi)


def _empty_plot(message: str, width: float, height: float, dpi: int, x_label: str = "", y_label: str = ""):
    """Create an empty plot with a centered message."""
    return (
        ggplot(pd.DataFrame({"x": [0], "y": [0]}), aes("x", "y"))
        + geom_text(aes(label=[message]), size=8)
        + labs(x=x_label, y=y_label)
        + _base_theme(width, height, dpi)
    )


def manhattan_plot(
    input_path: str,
    output_path: str,
    width: float = NPG_DOUBLE,
    height: float = 3.5,
    dpi: int = 300,
    log_scale: bool = False,
    label_top: int | None = None,
    highlight_genes: list[str] | None = None,
    max_ratio: float = 2.0,
    exclude_zero_ratio: bool = False,
    min_codons: int | None = None,
    min_variants: int | None = None,
) -> None:
    """Create Manhattan plot of per-gene piN/piS."""
    df = _load_gene_results(input_path)
    df = _apply_base_qc(df, min_codons=min_codons, min_variants=min_variants)
    df = _filter_ratio(df, max_ratio=max_ratio, exclude_zero_ratio=exclude_zero_ratio)

    if df.empty:
        _save_plot(_empty_plot("No genes with finite piN/piS", width, height, dpi, "Genomic position", "piN/piS"), output_path, dpi)
        return

    df, chrom_centers = _add_cumulative_pos(df)

    y_col = "piN_piS"
    y_label = "piN/piS"
    if log_scale:
        df["log2_piN_piS"] = np.log2(df["piN_piS"])
        y_col = "log2_piN_piS"
        y_label = "log\u2082(piN/piS)"

    neutral_y = 0.0 if log_scale else 1.0
    chroms = list(chrom_centers.keys())
    color_map = {c: _OKABE_ITO[i % len(_OKABE_ITO)] for i, c in enumerate(chroms)}

    p = (
        ggplot(df, aes(x="cum_pos", y=y_col, color="chrom"))
        + geom_point(size=0.8, alpha=0.7)
        + geom_hline(yintercept=neutral_y, linetype="dashed", color="red", alpha=0.6, size=0.3)
        + scale_color_manual(values=color_map)
        + guides(color=False)
        + scale_x_continuous(
            breaks=list(chrom_centers.values()),
            labels=chroms,
        )
        + labs(x="Chromosome", y=y_label, title="Per-gene piN/piS")
        + _base_theme(width, height, dpi)
        + theme(axis_text_x=element_text(rotation=45, ha="right", size=6))
    )

    label_df_parts = []
    if label_top is not None and label_top > 0:
        sorted_df = df.sort_values(y_col, ascending=False)
        label_df_parts.append(sorted_df.head(label_top))
    if highlight_genes:
        label_df_parts.append(df[df["gene_id"].isin(highlight_genes)])

    if label_df_parts:
        label_df = pd.concat(label_df_parts).drop_duplicates(subset=["gene_id"])
        p = p + geom_text(
            data=label_df,
            mapping=aes(x="cum_pos", y=y_col, label="gene_id"),
            size=5, nudge_y=0.05, color="black", inherit_aes=False,
        )

    _save_plot(p, output_path, dpi)


def scatter_plot(
    input_path: str,
    output_path: str,
    width: float = NPG_ONE_HALF,
    height: float = NPG_ONE_HALF,
    dpi: int = 300,
    color_by_chrom: bool = False,
    max_piN: float = 2.0,
    max_piS: float = 2.0,
    min_codons: int | None = None,
    min_variants: int | None = None,
) -> None:
    """Create piN vs piS scatter plot."""
    df = _load_gene_results(input_path)
    df = _apply_base_qc(df, min_codons=min_codons, min_variants=min_variants)
    df = _filter_metric(df, "piN", max_value=max_piN)
    df = _filter_metric(df, "piS", max_value=max_piS)

    if df.empty:
        _save_plot(_empty_plot("No genes with piN/piS data", width, height, dpi, "piS", "piN"), output_path, dpi)
        return

    chroms = _sort_chroms(df["chrom"].unique())
    df["chrom"] = pd.Categorical(df["chrom"], categories=chroms, ordered=True)

    mapping = aes(x="piS", y="piN", size="n_codons")
    if color_by_chrom:
        mapping = aes(x="piS", y="piN", size="n_codons", color="chrom")

    color_map = {c: _OKABE_ITO[i % len(_OKABE_ITO)] for i, c in enumerate(chroms)}

    p = (
        ggplot(df, mapping)
        + geom_point(alpha=0.6)
        + geom_abline(intercept=0, slope=1, linetype="dashed", color="red", alpha=0.6, size=0.3)
        + labs(x="piS", y="piN", title="piN vs piS", size="Codons")
        + _base_theme(width, height, dpi)
    )

    if color_by_chrom:
        p = p + scale_color_manual(values=color_map)

    _save_plot(p, output_path, dpi)


def histogram_plot(
    input_path: str,
    output_path: str,
    width: float = NPG_SINGLE,
    height: float = 3.0,
    dpi: int = 300,
    max_ratio: float = 2.0,
    exclude_zero_ratio: bool = False,
    min_codons: int | None = None,
    min_variants: int | None = None,
) -> None:
    """Create histogram of piN/piS distribution with density overlay."""
    df = _load_gene_results(input_path)
    df = _apply_base_qc(df, min_codons=min_codons, min_variants=min_variants)
    df = _filter_ratio(df, max_ratio=max_ratio, exclude_zero_ratio=exclude_zero_ratio)

    if df.empty:
        _save_plot(_empty_plot("No genes with finite piN/piS", width, height, dpi, "piN/piS", "Count"), output_path, dpi)
        return

    p = (
        ggplot(df, aes(x="piN_piS"))
        + geom_histogram(aes(y=after_stat("density")), bins=30, fill="#56B4E9", alpha=0.7, color="white")
        + geom_density(color="#0072B2", size=0.5)
        + geom_vline(xintercept=1.0, linetype="dashed", color="red", alpha=0.6, size=0.3)
        + labs(x="piN/piS", y="Density", title="Distribution of piN/piS")
        + _base_theme(width, height, dpi)
    )

    _save_plot(p, output_path, dpi)


def boxplot_plot(
    input_path: str,
    output_path: str,
    width: float = NPG_SINGLE,
    height: float = 7.0,
    dpi: int = 300,
    max_piN: float = 2.0,
    max_piS: float = 2.0,
    max_ratio: float = 2.0,
    exclude_zero_ratio: bool = False,
    min_codons: int | None = None,
    min_variants: int | None = None,
) -> None:
    """Create faceted boxplots of piN, piS, and piN/piS per chromosome."""
    df = _load_gene_results(input_path)
    df = _apply_base_qc(df, min_codons=min_codons, min_variants=min_variants)
    chroms = _sort_chroms(df["chrom"].unique())
    df["chrom"] = pd.Categorical(df["chrom"], categories=chroms, ordered=True)

    # Per-facet filtering: apply metric-specific filters independently
    df_piN = _filter_metric(df, "piN", max_value=max_piN)[["chrom", "gene_id", "piN"]].rename(columns={"piN": "value"})
    df_piN["metric"] = "piN"
    df_piS = _filter_metric(df, "piS", max_value=max_piS)[["chrom", "gene_id", "piS"]].rename(columns={"piS": "value"})
    df_piS["metric"] = "piS"
    df_ratio = _filter_ratio(df, max_ratio=max_ratio, exclude_zero_ratio=exclude_zero_ratio)[["chrom", "gene_id", "piN_piS"]].rename(columns={"piN_piS": "value"})
    df_ratio["metric"] = "piN_piS"
    long = pd.concat([df_piN, df_piS, df_ratio], ignore_index=True)

    if long.empty:
        _save_plot(_empty_plot("No data to plot", width, height, dpi), output_path, dpi)
        return

    long["metric"] = pd.Categorical(
        long["metric"],
        categories=["piN", "piS", "piN_piS"],
        ordered=True,
    )

    color_map = {c: _OKABE_ITO[i % len(_OKABE_ITO)] for i, c in enumerate(chroms)}

    p = (
        ggplot(long, aes(x="chrom", y="value", fill="chrom"))
        + geom_boxplot(alpha=0.7, show_legend=False)
        + facet_wrap("metric", scales="free_y", ncol=1)
        + scale_fill_manual(values=color_map)
        + labs(x="Chromosome", y="", title="piN, piS, piN/piS by chromosome")
        + _base_theme(width, height, dpi)
        + theme(axis_text_x=element_text(rotation=45, ha="right", size=6))
    )

    _save_plot(p, output_path, dpi)


def sliding_window_plot(
    input_path: str,
    output_path: str,
    width: float = NPG_DOUBLE,
    height: float | None = None,
    dpi: int = 300,
    max_ratio: float = 2.0,
) -> None:
    """Create sliding window piN/piS line plot faceted by chromosome.

    Height defaults to 1.5 inches per chromosome when not specified.
    """
    df = pd.read_csv(input_path, sep="\t")
    df["chrom"] = df["chrom"].astype(str)
    df = df.dropna(subset=["piN_piS"])
    df = df[(df["piN_piS"] >= 0) & (df["piN_piS"] <= max_ratio)]

    if df.empty:
        h = height or 4.0
        _save_plot(_empty_plot("No sliding window data", width, h, dpi, "Position", "piN/piS"), output_path, dpi)
        return

    chroms = _sort_chroms(df["chrom"].unique())
    df["chrom"] = pd.Categorical(df["chrom"], categories=chroms, ordered=True)
    df["win_mid"] = (df["win_start"] + df["win_end"]) / 2

    if height is None:
        height = 1.5 * len(chroms) + 0.5  # 1.5 in per panel + padding

    p = (
        ggplot(df, aes(x="win_mid", y="piN_piS", group="gene_id"))
        + geom_line(color="#0072B2", alpha=0.8, size=0.3)
        + geom_hline(yintercept=1.0, linetype="dashed", color="red", alpha=0.6, size=0.3)
        + facet_wrap("chrom", scales="free_x", ncol=1)
        + labs(x="Genomic position (bp)", y="piN/piS", title="Sliding window piN/piS")
        + _base_theme(width, height, dpi)
    )

    _save_plot(p, output_path, dpi)
