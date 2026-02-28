"""Manhattan plot for per-gene piN/piS."""

import matplotlib
matplotlib.use("Agg")  # non-interactive backend
import matplotlib.pyplot as plt
import pandas as pd


def manhattan_plot(
    gene_results_path: str,
    output_path: str,
    width: float = 16.0,
    height: float = 6.0,
) -> None:
    """Create Manhattan plot of per-gene piN/piS.

    Args:
        gene_results_path: Path to gene_results.tsv
        output_path: Path for output PNG
        width: Figure width in inches
        height: Figure height in inches
    """
    df = pd.read_csv(gene_results_path, sep="\t")

    # Drop genes with NA piN_piS
    df = df.dropna(subset=["piN_piS"])
    if df.empty:
        # Create empty plot with message
        fig, ax = plt.subplots(figsize=(width, height))
        ax.text(
            0.5, 0.5, "No genes with finite piN/piS",
            transform=ax.transAxes, ha="center", va="center", fontsize=14,
        )
        ax.set_xlabel("Genomic position")
        ax.set_ylabel("piN/piS")
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return

    # Calculate cumulative positions per chromosome
    chroms = df["chrom"].unique()
    # Sort chromosomes: numeric first, then alpha
    chroms = sorted(
        chroms,
        key=lambda c: (
            not str(c).replace("chr", "").isdigit(),
            int(str(c).replace("chr", "")) if str(c).replace("chr", "").isdigit() else str(c),
        ),
    )

    chrom_offsets = {}
    cumulative = 0
    chrom_centers = {}
    for chrom in chroms:
        chrom_offsets[chrom] = cumulative
        chrom_data = df[df["chrom"] == chrom]
        chrom_max = chrom_data["end"].max()
        chrom_centers[chrom] = cumulative + chrom_max / 2
        cumulative += chrom_max + chrom_max * 0.05  # 5% gap

    # Gene midpoint in cumulative coordinates
    df["cum_pos"] = df.apply(
        lambda row: chrom_offsets[row["chrom"]] + (row["start"] + row["end"]) / 2,
        axis=1,
    )

    # Alternate colors by chromosome
    colors = ["#1f77b4", "#aec7e8"]  # dark blue, light blue

    fig, ax = plt.subplots(figsize=(width, height))

    for i, chrom in enumerate(chroms):
        mask = df["chrom"] == chrom
        ax.scatter(
            df.loc[mask, "cum_pos"],
            df.loc[mask, "piN_piS"],
            c=colors[i % 2],
            s=8,
            alpha=0.7,
            label=chrom,
            edgecolors="none",
        )

    # Neutral expectation line
    ax.axhline(y=1.0, color="red", linestyle="--", linewidth=0.8, alpha=0.6)

    # Chromosome labels
    ax.set_xticks([chrom_centers[c] for c in chroms])
    ax.set_xticklabels(chroms, rotation=45, ha="right", fontsize=8)

    ax.set_xlabel("Chromosome")
    ax.set_ylabel("piN/piS")
    ax.set_title("Per-gene piN/piS Manhattan Plot")

    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
