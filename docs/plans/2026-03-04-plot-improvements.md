# Plot Improvements Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Replace the minimal matplotlib Manhattan plot with a full suite of 5 publication-ready plot types using plotnine, exposed as `pie plot <subcommand>`.

**Architecture:** Single `plot.py` module with plotnine. Shared helpers for theming, data loading, and chromosome sorting. CLI restructured so `pie plot` is a `click.Group` with subcommands: `manhattan`, `scatter`, `histogram`, `boxplot`, `sliding-window`. Each subcommand has shared options (`-i`, `-o`, `-W`, `-H`, `--dpi`) plus type-specific flags.

**Tech Stack:** plotnine (ggplot2 for Python), pandas, click

**Design doc:** `docs/plans/2026-03-04-plot-improvements-design.md`

---

### Task 1: Add plotnine dependency

**Files:**
- Modify: `environment.yml`
- Modify: `pyproject.toml`

**Step 1: Update environment.yml**

Replace `matplotlib` with `plotnine` (plotnine depends on matplotlib, so it's still available transitively):

```yaml
# environment.yml
name: pie
channels:
  - conda-forge
  - bioconda
dependencies:
  - python=3.12
  - cyvcf2
  - gffutils
  - pysam
  - numpy
  - pandas
  - plotnine
  - click
  - pytest
  - htslib
```

**Step 2: Update pyproject.toml**

Replace `"matplotlib"` with `"plotnine"` in the dependencies list:

```toml
dependencies = [
    "cyvcf2",
    "gffutils",
    "pysam",
    "numpy",
    "pandas",
    "plotnine",
    "click",
]
```

**Step 3: Install plotnine**

Run: `pip install plotnine`

**Step 4: Verify import works**

Run: `python -c "from plotnine import ggplot, aes, geom_point; print('ok')"`
Expected: `ok`

**Step 5: Commit**

```bash
git add environment.yml pyproject.toml
git commit -m "Replace matplotlib with plotnine dependency"
```

---

### Task 2: Rewrite plot.py with shared helpers and manhattan_plot

**Files:**
- Modify: `src/pie/plot.py`
- Modify: `tests/test_plot.py`

**Step 1: Write the failing tests**

Replace `tests/test_plot.py` with tests for the new plotnine-based manhattan_plot. The function signature changes: it now accepts `log_scale`, `label_top`, and `highlight_genes` parameters, and saves via plotnine's `ggsave`.

```python
import pandas as pd
import pytest
from pie.plot import manhattan_plot


def _make_gene_tsv(tmp_path, df):
    """Helper: write a DataFrame to gene_results.tsv and return the path string."""
    tsv = tmp_path / "gene_results.tsv"
    df.to_csv(tsv, sep="\t", index=False)
    return str(tsv)


@pytest.fixture
def basic_gene_df():
    return pd.DataFrame({
        "chrom": ["chr1", "chr1", "chr2", "chr2"],
        "gene_id": ["g1", "g2", "g3", "g4"],
        "start": [100, 500, 100, 800],
        "end": [200, 600, 300, 900],
        "n_codons": [30, 40, 25, 50],
        "piN": [0.01, 0.02, 0.015, 0.03],
        "piS": [0.02, 0.017, 0.019, 0.015],
        "piN_piS": [0.5, 1.2, 0.8, 2.0],
    })


class TestManhattanPlot:
    def test_creates_png(self, tmp_path, basic_gene_df):
        tsv = _make_gene_tsv(tmp_path, basic_gene_df)
        out = str(tmp_path / "manhattan.png")
        manhattan_plot(tsv, out)
        assert (tmp_path / "manhattan.png").exists()
        assert (tmp_path / "manhattan.png").stat().st_size > 0

    def test_custom_dimensions(self, tmp_path, basic_gene_df):
        tsv = _make_gene_tsv(tmp_path, basic_gene_df)
        out = str(tmp_path / "manhattan.png")
        manhattan_plot(tsv, out, width=20, height=8)
        assert (tmp_path / "manhattan.png").exists()

    def test_log_scale(self, tmp_path, basic_gene_df):
        tsv = _make_gene_tsv(tmp_path, basic_gene_df)
        out = str(tmp_path / "manhattan.png")
        manhattan_plot(tsv, out, log_scale=True)
        assert (tmp_path / "manhattan.png").exists()

    def test_label_top(self, tmp_path, basic_gene_df):
        tsv = _make_gene_tsv(tmp_path, basic_gene_df)
        out = str(tmp_path / "manhattan.png")
        manhattan_plot(tsv, out, label_top=2)
        assert (tmp_path / "manhattan.png").exists()

    def test_highlight_genes(self, tmp_path, basic_gene_df):
        tsv = _make_gene_tsv(tmp_path, basic_gene_df)
        out = str(tmp_path / "manhattan.png")
        manhattan_plot(tsv, out, highlight_genes=["g1", "g3"])
        assert (tmp_path / "manhattan.png").exists()

    def test_numeric_chrom_ids(self, tmp_path):
        df = pd.DataFrame({
            "chrom": [1, 1, 2, 2],
            "gene_id": ["g1", "g2", "g3", "g4"],
            "start": [100, 500, 100, 800],
            "end": [200, 600, 300, 900],
            "n_codons": [30, 40, 25, 50],
            "piN": [0.01, 0.02, 0.015, 0.03],
            "piS": [0.02, 0.017, 0.019, 0.015],
            "piN_piS": [0.5, 1.2, 0.8, 2.0],
        })
        tsv = _make_gene_tsv(tmp_path, df)
        out = str(tmp_path / "manhattan.png")
        manhattan_plot(tsv, out)
        assert (tmp_path / "manhattan.png").exists()

    def test_handles_all_na(self, tmp_path):
        df = pd.DataFrame({
            "chrom": ["chr1"],
            "gene_id": ["g1"],
            "start": [100],
            "end": [200],
            "n_codons": [30],
            "piN": [0.0],
            "piS": [0.0],
            "piN_piS": [float("nan")],
        })
        tsv = _make_gene_tsv(tmp_path, df)
        out = str(tmp_path / "manhattan.png")
        manhattan_plot(tsv, out)
        assert (tmp_path / "manhattan.png").exists()

    def test_pdf_output(self, tmp_path, basic_gene_df):
        tsv = _make_gene_tsv(tmp_path, basic_gene_df)
        out = str(tmp_path / "manhattan.pdf")
        manhattan_plot(tsv, out)
        assert (tmp_path / "manhattan.pdf").exists()
```

**Step 2: Run tests to verify they fail**

Run: `pytest tests/test_plot.py -v`
Expected: FAIL (old manhattan_plot signature doesn't accept new params)

**Step 3: Write the implementation**

Replace `src/pie/plot.py` entirely:

```python
"""Publication-ready plots for per-gene and sliding-window piN/piS results."""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
from plotnine import (
    aes,
    element_blank,
    element_text,
    geom_hline,
    geom_point,
    geom_text,
    ggplot,
    ggsave,
    labs,
    scale_color_manual,
    scale_x_continuous,
    theme,
    theme_bw,
)

# Okabe-Ito colorblind-safe palette (8 colors, cycles for >8 chroms)
_OKABE_ITO = [
    "#E69F00", "#56B4E9", "#009E73", "#F0E442",
    "#0072B2", "#D55E00", "#CC79A7", "#999999",
]


def _base_theme(width: float = 12, height: float = 6, dpi: int = 300) -> theme:
    """Publication-quality base theme."""
    return (
        theme_bw()
        + theme(
            figure_size=(width, height),
            dpi=dpi,
            axis_text=element_text(size=9),
            axis_title=element_text(size=11),
            plot_title=element_text(size=13, weight="bold"),
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


def _load_gene_results(path: str) -> pd.DataFrame:
    """Load gene_results.tsv and prepare for plotting."""
    df = pd.read_csv(path, sep="\t")
    df["chrom"] = df["chrom"].astype(str)
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
        cumulative += chrom_max * 1.05  # 5% gap

    df = df.copy()
    df["cum_pos"] = df.apply(
        lambda row: chrom_offsets[row["chrom"]] + (row["start"] + row["end"]) / 2,
        axis=1,
    )
    # Assign color index by chromosome order
    chrom_to_idx = {c: i for i, c in enumerate(chroms)}
    df["chrom_color"] = df["chrom"].map(lambda c: _OKABE_ITO[chrom_to_idx[c] % len(_OKABE_ITO)])
    df["chrom"] = pd.Categorical(df["chrom"], categories=chroms, ordered=True)
    return df, chrom_centers


def manhattan_plot(
    input_path: str,
    output_path: str,
    width: float = 12,
    height: float = 6,
    dpi: int = 300,
    log_scale: bool = False,
    label_top: int | None = None,
    highlight_genes: list[str] | None = None,
) -> None:
    """Create Manhattan plot of per-gene piN/piS."""
    df = _load_gene_results(input_path)
    df = df.dropna(subset=["piN_piS"])

    if df.empty:
        # Create minimal empty plot
        p = (
            ggplot(pd.DataFrame({"x": [0], "y": [0]}), aes("x", "y"))
            + geom_text(aes(label=["No genes with finite piN/piS"]), size=14)
            + labs(x="Genomic position", y="piN/piS")
            + _base_theme(width, height, dpi)
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ggsave(p, output_path, dpi=dpi)
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
        + geom_point(size=1.5, alpha=0.7)
        + geom_hline(yintercept=neutral_y, linetype="dashed", color="red", alpha=0.6, size=0.5)
        + scale_color_manual(values=color_map, guide=False)
        + scale_x_continuous(
            breaks=list(chrom_centers.values()),
            labels=chroms,
        )
        + labs(x="Chromosome", y=y_label, title="Per-gene piN/piS")
        + _base_theme(width, height, dpi)
        + theme(axis_text_x=element_text(rotation=45, ha="right", size=8))
    )

    # Label outlier genes
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
            size=7, nudge_y=0.1, color="black", inherit_aes=False,
        )

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ggsave(p, output_path, dpi=dpi)
```

**Step 4: Run tests to verify they pass**

Run: `pytest tests/test_plot.py -v`
Expected: All PASS

**Step 5: Commit**

```bash
git add src/pie/plot.py tests/test_plot.py
git commit -m "Rewrite manhattan_plot with plotnine and shared helpers"
```

---

### Task 3: Add scatter_plot

**Files:**
- Modify: `src/pie/plot.py`
- Modify: `tests/test_plot.py`

**Step 1: Write the failing tests**

Append to `tests/test_plot.py`:

```python
from pie.plot import scatter_plot


class TestScatterPlot:
    def test_creates_png(self, tmp_path, basic_gene_df):
        tsv = _make_gene_tsv(tmp_path, basic_gene_df)
        out = str(tmp_path / "scatter.png")
        scatter_plot(tsv, out)
        assert (tmp_path / "scatter.png").exists()
        assert (tmp_path / "scatter.png").stat().st_size > 0

    def test_color_by_chrom(self, tmp_path, basic_gene_df):
        tsv = _make_gene_tsv(tmp_path, basic_gene_df)
        out = str(tmp_path / "scatter.png")
        scatter_plot(tsv, out, color_by_chrom=True)
        assert (tmp_path / "scatter.png").exists()

    def test_handles_all_na(self, tmp_path):
        df = pd.DataFrame({
            "chrom": ["chr1"],
            "gene_id": ["g1"],
            "start": [100],
            "end": [200],
            "n_codons": [30],
            "piN": [float("nan")],
            "piS": [float("nan")],
            "piN_piS": [float("nan")],
        })
        tsv = _make_gene_tsv(tmp_path, df)
        out = str(tmp_path / "scatter.png")
        scatter_plot(tsv, out)
        assert (tmp_path / "scatter.png").exists()
```

**Step 2: Run tests to verify they fail**

Run: `pytest tests/test_plot.py::TestScatterPlot -v`
Expected: FAIL (scatter_plot not defined)

**Step 3: Write minimal implementation**

Add to `src/pie/plot.py`:

```python
from plotnine import geom_abline

def scatter_plot(
    input_path: str,
    output_path: str,
    width: float = 12,
    height: float = 6,
    dpi: int = 300,
    color_by_chrom: bool = False,
) -> None:
    """Create piN vs piS scatter plot."""
    df = _load_gene_results(input_path)
    df = df.dropna(subset=["piN", "piS"])

    if df.empty:
        p = (
            ggplot(pd.DataFrame({"x": [0], "y": [0]}), aes("x", "y"))
            + geom_text(aes(label=["No genes with piN/piS data"]), size=14)
            + labs(x="piS", y="piN")
            + _base_theme(width, height, dpi)
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ggsave(p, output_path, dpi=dpi)
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
        + geom_abline(intercept=0, slope=1, linetype="dashed", color="red", alpha=0.6, size=0.5)
        + labs(x="piS", y="piN", title="piN vs piS", size="Codons")
        + _base_theme(width, height, dpi)
    )

    if color_by_chrom:
        p = p + scale_color_manual(values=color_map)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ggsave(p, output_path, dpi=dpi)
```

Add `geom_abline` to the imports at the top of `plot.py`.

**Step 4: Run tests to verify they pass**

Run: `pytest tests/test_plot.py::TestScatterPlot -v`
Expected: All PASS

**Step 5: Commit**

```bash
git add src/pie/plot.py tests/test_plot.py
git commit -m "Add piN vs piS scatter plot"
```

---

### Task 4: Add histogram_plot

**Files:**
- Modify: `src/pie/plot.py`
- Modify: `tests/test_plot.py`

**Step 1: Write the failing tests**

Append to `tests/test_plot.py`:

```python
from pie.plot import histogram_plot


class TestHistogramPlot:
    def test_creates_png(self, tmp_path, basic_gene_df):
        tsv = _make_gene_tsv(tmp_path, basic_gene_df)
        out = str(tmp_path / "histogram.png")
        histogram_plot(tsv, out)
        assert (tmp_path / "histogram.png").exists()
        assert (tmp_path / "histogram.png").stat().st_size > 0

    def test_handles_all_na(self, tmp_path):
        df = pd.DataFrame({
            "chrom": ["chr1"],
            "gene_id": ["g1"],
            "start": [100],
            "end": [200],
            "n_codons": [30],
            "piN": [0.0],
            "piS": [0.0],
            "piN_piS": [float("nan")],
        })
        tsv = _make_gene_tsv(tmp_path, df)
        out = str(tmp_path / "histogram.png")
        histogram_plot(tsv, out)
        assert (tmp_path / "histogram.png").exists()
```

**Step 2: Run tests to verify they fail**

Run: `pytest tests/test_plot.py::TestHistogramPlot -v`
Expected: FAIL (histogram_plot not defined)

**Step 3: Write minimal implementation**

Add to `src/pie/plot.py`:

```python
from plotnine import geom_histogram, geom_density, geom_vline

def histogram_plot(
    input_path: str,
    output_path: str,
    width: float = 12,
    height: float = 6,
    dpi: int = 300,
) -> None:
    """Create histogram of piN/piS distribution with density overlay."""
    df = _load_gene_results(input_path)
    df = df.dropna(subset=["piN_piS"])

    if df.empty:
        p = (
            ggplot(pd.DataFrame({"x": [0], "y": [0]}), aes("x", "y"))
            + geom_text(aes(label=["No genes with finite piN/piS"]), size=14)
            + labs(x="piN/piS", y="Count")
            + _base_theme(width, height, dpi)
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ggsave(p, output_path, dpi=dpi)
        return

    p = (
        ggplot(df, aes(x="piN_piS"))
        + geom_histogram(aes(y="..density.."), bins=30, fill="#56B4E9", alpha=0.7, color="white")
        + geom_density(color="#0072B2", size=1)
        + geom_vline(xintercept=1.0, linetype="dashed", color="red", alpha=0.6, size=0.5)
        + labs(x="piN/piS", y="Density", title="Distribution of piN/piS")
        + _base_theme(width, height, dpi)
    )

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ggsave(p, output_path, dpi=dpi)
```

Add `geom_histogram`, `geom_density`, `geom_vline` to the imports.

**Step 4: Run tests to verify they pass**

Run: `pytest tests/test_plot.py::TestHistogramPlot -v`
Expected: All PASS

**Step 5: Commit**

```bash
git add src/pie/plot.py tests/test_plot.py
git commit -m "Add piN/piS histogram with density overlay"
```

---

### Task 5: Add boxplot_plot

**Files:**
- Modify: `src/pie/plot.py`
- Modify: `tests/test_plot.py`

**Step 1: Write the failing tests**

Append to `tests/test_plot.py`:

```python
from pie.plot import boxplot_plot


class TestBoxplotPlot:
    def test_creates_png(self, tmp_path, basic_gene_df):
        tsv = _make_gene_tsv(tmp_path, basic_gene_df)
        out = str(tmp_path / "boxplot.png")
        boxplot_plot(tsv, out)
        assert (tmp_path / "boxplot.png").exists()
        assert (tmp_path / "boxplot.png").stat().st_size > 0

    def test_single_chromosome(self, tmp_path):
        df = pd.DataFrame({
            "chrom": ["chr1", "chr1", "chr1"],
            "gene_id": ["g1", "g2", "g3"],
            "start": [100, 300, 500],
            "end": [200, 400, 600],
            "n_codons": [30, 40, 25],
            "piN": [0.01, 0.02, 0.015],
            "piS": [0.02, 0.017, 0.019],
            "piN_piS": [0.5, 1.2, 0.8],
        })
        tsv = _make_gene_tsv(tmp_path, df)
        out = str(tmp_path / "boxplot.png")
        boxplot_plot(tsv, out)
        assert (tmp_path / "boxplot.png").exists()

    def test_handles_all_na(self, tmp_path):
        df = pd.DataFrame({
            "chrom": ["chr1"],
            "gene_id": ["g1"],
            "start": [100],
            "end": [200],
            "n_codons": [30],
            "piN": [0.0],
            "piS": [0.0],
            "piN_piS": [float("nan")],
        })
        tsv = _make_gene_tsv(tmp_path, df)
        out = str(tmp_path / "boxplot.png")
        boxplot_plot(tsv, out)
        assert (tmp_path / "boxplot.png").exists()
```

**Step 2: Run tests to verify they fail**

Run: `pytest tests/test_plot.py::TestBoxplotPlot -v`
Expected: FAIL (boxplot_plot not defined)

**Step 3: Write minimal implementation**

Add to `src/pie/plot.py`:

```python
from plotnine import facet_wrap, geom_boxplot

def boxplot_plot(
    input_path: str,
    output_path: str,
    width: float = 12,
    height: float = 8,
    dpi: int = 300,
) -> None:
    """Create faceted boxplots of piN, piS, and piN/piS per chromosome."""
    df = _load_gene_results(input_path)
    chroms = _sort_chroms(df["chrom"].unique())
    df["chrom"] = pd.Categorical(df["chrom"], categories=chroms, ordered=True)

    # Melt piN, piS, piN_piS into long format for faceting
    long = df.melt(
        id_vars=["chrom", "gene_id"],
        value_vars=["piN", "piS", "piN_piS"],
        var_name="metric",
        value_name="value",
    ).dropna(subset=["value"])

    if long.empty:
        p = (
            ggplot(pd.DataFrame({"x": [0], "y": [0]}), aes("x", "y"))
            + geom_text(aes(label=["No data to plot"]), size=14)
            + _base_theme(width, height, dpi)
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ggsave(p, output_path, dpi=dpi)
        return

    # Order metrics for display
    long["metric"] = pd.Categorical(
        long["metric"],
        categories=["piN", "piS", "piN_piS"],
        ordered=True,
    )

    p = (
        ggplot(long, aes(x="chrom", y="value", fill="chrom"))
        + geom_boxplot(alpha=0.7, show_legend=False)
        + facet_wrap("metric", scales="free_y", ncol=1)
        + labs(x="Chromosome", y="", title="piN, piS, piN/piS by chromosome")
        + _base_theme(width, height, dpi)
        + theme(axis_text_x=element_text(rotation=45, ha="right", size=8))
    )

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ggsave(p, output_path, dpi=dpi)
```

Add `facet_wrap`, `geom_boxplot`, and `scale_fill_manual` to the imports. Also use `scale_fill_manual` with the Okabe-Ito palette:

```python
    color_map = {c: _OKABE_ITO[i % len(_OKABE_ITO)] for i, c in enumerate(chroms)}
    # Add after geom_boxplot line:
    + scale_fill_manual(values=color_map)
```

**Step 4: Run tests to verify they pass**

Run: `pytest tests/test_plot.py::TestBoxplotPlot -v`
Expected: All PASS

**Step 5: Commit**

```bash
git add src/pie/plot.py tests/test_plot.py
git commit -m "Add faceted boxplot of piN/piS/piN_piS per chromosome"
```

---

### Task 6: Add sliding_window_plot

**Files:**
- Modify: `src/pie/plot.py`
- Modify: `tests/test_plot.py`

**Step 1: Write the failing tests**

Append to `tests/test_plot.py`:

```python
from pie.plot import sliding_window_plot


class TestSlidingWindowPlot:
    def test_creates_png(self, tmp_path):
        df = pd.DataFrame({
            "chrom": ["chr1"] * 5 + ["chr2"] * 5,
            "win_start": list(range(0, 500, 100)) * 2,
            "win_end": list(range(100, 600, 100)) * 2,
            "gene_id": ["g1"] * 5 + ["g2"] * 5,
            "piN_piS": [0.5, 0.8, 1.0, 1.2, 0.9, 1.5, 1.1, 0.7, 0.6, 1.3],
        })
        tsv = tmp_path / "window_results.tsv"
        df.to_csv(tsv, sep="\t", index=False)
        out = str(tmp_path / "sliding_window.png")
        sliding_window_plot(str(tsv), out)
        assert (tmp_path / "sliding_window.png").exists()
        assert (tmp_path / "sliding_window.png").stat().st_size > 0

    def test_handles_empty(self, tmp_path):
        df = pd.DataFrame({
            "chrom": pd.Series([], dtype=str),
            "win_start": pd.Series([], dtype=int),
            "win_end": pd.Series([], dtype=int),
            "gene_id": pd.Series([], dtype=str),
            "piN_piS": pd.Series([], dtype=float),
        })
        tsv = tmp_path / "window_results.tsv"
        df.to_csv(tsv, sep="\t", index=False)
        out = str(tmp_path / "sliding_window.png")
        sliding_window_plot(str(tsv), out)
        assert (tmp_path / "sliding_window.png").exists()
```

**Step 2: Run tests to verify they fail**

Run: `pytest tests/test_plot.py::TestSlidingWindowPlot -v`
Expected: FAIL (sliding_window_plot not defined)

**Step 3: Write minimal implementation**

Add to `src/pie/plot.py`:

```python
from plotnine import geom_line

def sliding_window_plot(
    input_path: str,
    output_path: str,
    width: float = 12,
    height: float = 6,
    dpi: int = 300,
) -> None:
    """Create sliding window piN/piS line plot faceted by chromosome."""
    df = pd.read_csv(input_path, sep="\t")
    df["chrom"] = df["chrom"].astype(str)
    df = df.dropna(subset=["piN_piS"])

    if df.empty:
        p = (
            ggplot(pd.DataFrame({"x": [0], "y": [0]}), aes("x", "y"))
            + geom_text(aes(label=["No sliding window data"]), size=14)
            + labs(x="Position", y="piN/piS")
            + _base_theme(width, height, dpi)
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ggsave(p, output_path, dpi=dpi)
        return

    chroms = _sort_chroms(df["chrom"].unique())
    df["chrom"] = pd.Categorical(df["chrom"], categories=chroms, ordered=True)
    df["win_mid"] = (df["win_start"] + df["win_end"]) / 2

    p = (
        ggplot(df, aes(x="win_mid", y="piN_piS"))
        + geom_line(color="#0072B2", alpha=0.8, size=0.5)
        + geom_hline(yintercept=1.0, linetype="dashed", color="red", alpha=0.6, size=0.5)
        + facet_wrap("chrom", scales="free_x")
        + labs(x="Genomic position (bp)", y="piN/piS", title="Sliding window piN/piS")
        + _base_theme(width, height, dpi)
    )

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ggsave(p, output_path, dpi=dpi)
```

Add `geom_line` to the imports.

**Step 4: Run tests to verify they pass**

Run: `pytest tests/test_plot.py::TestSlidingWindowPlot -v`
Expected: All PASS

**Step 5: Commit**

```bash
git add src/pie/plot.py tests/test_plot.py
git commit -m "Add sliding window piN/piS line plot"
```

---

### Task 7: Restructure CLI — pie plot as click.Group with subcommands

**Files:**
- Modify: `src/pie/cli.py`
- Modify: `tests/test_cli.py`

**Step 1: Write the failing tests**

Replace `TestPlotCommand` in `tests/test_cli.py` with:

```python
class TestPlotCommand:
    def test_plot_group_help(self, runner):
        result = runner.invoke(main, ["plot", "--help"])
        assert result.exit_code == 0
        assert "manhattan" in result.output
        assert "scatter" in result.output
        assert "histogram" in result.output
        assert "boxplot" in result.output
        assert "sliding-window" in result.output

    def test_manhattan_help(self, runner):
        result = runner.invoke(main, ["plot", "manhattan", "--help"])
        assert result.exit_code == 0
        assert "--input" in result.output
        assert "--output" in result.output
        assert "--log-scale" in result.output
        assert "--label-top" in result.output
        assert "--highlight-genes" in result.output

    def test_scatter_help(self, runner):
        result = runner.invoke(main, ["plot", "scatter", "--help"])
        assert result.exit_code == 0
        assert "--color-by-chrom" in result.output

    def test_manhattan_from_results(self, runner, ref_fasta, gff3_file, vcf_file, tmp_path):
        runner.invoke(main, [
            "run", "--vcf", vcf_file, "--gff", gff3_file, "--fasta", ref_fasta,
            "--outdir", str(tmp_path), "--min-freq", "0", "--min-depth", "0", "--min-qual", "0",
        ])
        result = runner.invoke(main, [
            "plot", "manhattan",
            "-i", str(tmp_path / "gene_results.tsv"),
            "-o", str(tmp_path / "manhattan.png"),
        ])
        assert result.exit_code == 0, result.output
        assert (tmp_path / "manhattan.png").exists()

    def test_scatter_from_results(self, runner, ref_fasta, gff3_file, vcf_file, tmp_path):
        runner.invoke(main, [
            "run", "--vcf", vcf_file, "--gff", gff3_file, "--fasta", ref_fasta,
            "--outdir", str(tmp_path), "--min-freq", "0", "--min-depth", "0", "--min-qual", "0",
        ])
        result = runner.invoke(main, [
            "plot", "scatter",
            "-i", str(tmp_path / "gene_results.tsv"),
            "-o", str(tmp_path / "scatter.png"),
        ])
        assert result.exit_code == 0, result.output
        assert (tmp_path / "scatter.png").exists()

    def test_boxplot_from_results(self, runner, ref_fasta, gff3_file, vcf_file, tmp_path):
        runner.invoke(main, [
            "run", "--vcf", vcf_file, "--gff", gff3_file, "--fasta", ref_fasta,
            "--outdir", str(tmp_path), "--min-freq", "0", "--min-depth", "0", "--min-qual", "0",
        ])
        result = runner.invoke(main, [
            "plot", "boxplot",
            "-i", str(tmp_path / "gene_results.tsv"),
            "-o", str(tmp_path / "boxplot.png"),
        ])
        assert result.exit_code == 0, result.output
        assert (tmp_path / "boxplot.png").exists()

    def test_histogram_from_results(self, runner, ref_fasta, gff3_file, vcf_file, tmp_path):
        runner.invoke(main, [
            "run", "--vcf", vcf_file, "--gff", gff3_file, "--fasta", ref_fasta,
            "--outdir", str(tmp_path), "--min-freq", "0", "--min-depth", "0", "--min-qual", "0",
        ])
        result = runner.invoke(main, [
            "plot", "histogram",
            "-i", str(tmp_path / "gene_results.tsv"),
            "-o", str(tmp_path / "histogram.png"),
        ])
        assert result.exit_code == 0, result.output
        assert (tmp_path / "histogram.png").exists()

    def test_sliding_window_from_results(self, runner, ref_fasta, gff3_file, vcf_file, tmp_path):
        runner.invoke(main, [
            "run", "--vcf", vcf_file, "--gff", gff3_file, "--fasta", ref_fasta,
            "--outdir", str(tmp_path), "--min-freq", "0", "--min-depth", "0", "--min-qual", "0",
        ])
        result = runner.invoke(main, [
            "plot", "sliding-window",
            "-i", str(tmp_path / "window_results.tsv"),
            "-o", str(tmp_path / "sw.png"),
        ])
        assert result.exit_code == 0, result.output
        assert (tmp_path / "sw.png").exists()
```

**Step 2: Run tests to verify they fail**

Run: `pytest tests/test_cli.py::TestPlotCommand -v`
Expected: FAIL (old `pie plot` is not a group, doesn't accept subcommands)

**Step 3: Write the implementation**

Replace the `plot` command in `src/pie/cli.py` (lines 293-318) with:

```python
@main.group(invoke_without_command=True, no_args_is_help=True, context_settings=_HELP_OPTS)
def plot():
    """Create publication-ready plots from piN/piS results.

    \b
    Subcommands:
      manhattan       Genome-wide Manhattan plot of per-gene piN/piS
      scatter         piN vs piS scatter plot
      histogram       Distribution of piN/piS ratios
      boxplot         Per-chromosome boxplots of piN, piS, piN/piS
      sliding-window  Sliding window piN/piS along chromosomes
    """


def _shared_plot_options(f):
    """Decorator: shared options for all plot subcommands."""
    f = click.option("-i", "--input", "input_path", required=True, help="Input TSV file.")(f)
    f = click.option("-o", "--output", "output_path", required=True, help="Output plot path (PNG/PDF/SVG).")(f)
    f = click.option("-W", "--width", default=12.0, show_default=True, help="Figure width in inches.")(f)
    f = click.option("-H", "--height", default=6.0, show_default=True, help="Figure height in inches.")(f)
    f = click.option("--dpi", default=300, show_default=True, help="Resolution in dots per inch.")(f)
    return f


@plot.command(no_args_is_help=True, context_settings=_HELP_OPTS)
@_shared_plot_options
@click.option("--log-scale", is_flag=True, help="Use log2 scale for piN/piS y-axis.")
@click.option("--label-top", default=None, type=int, help="Label top N outlier genes.")
@click.option("--highlight-genes", default=None, help="Comma-separated gene IDs to label.")
def manhattan(input_path, output_path, width, height, dpi, log_scale, label_top, highlight_genes):
    """Create Manhattan plot of per-gene piN/piS.

    \b
    Genes are ordered by chromosomal position along the x-axis.  A dashed
    line marks neutral expectation (piN/piS = 1).

    \b
    Examples:
      pie plot manhattan -i gene_results.tsv -o manhattan.png
      pie plot manhattan -i gene_results.tsv -o manhattan.png --log-scale
      pie plot manhattan -i gene_results.tsv -o manhattan.png --label-top 10
    """
    from pie.plot import manhattan_plot as _manhattan_plot

    if not os.path.exists(input_path):
        click.echo(f"Error: file not found: {input_path}", err=True)
        sys.exit(1)

    genes = [g.strip() for g in highlight_genes.split(",")] if highlight_genes else None
    _manhattan_plot(input_path, output_path, width=width, height=height, dpi=dpi,
                    log_scale=log_scale, label_top=label_top, highlight_genes=genes)
    click.echo(f"Plot saved to {output_path}")


@plot.command(no_args_is_help=True, context_settings=_HELP_OPTS)
@_shared_plot_options
@click.option("--color-by-chrom", is_flag=True, help="Color points by chromosome.")
def scatter(input_path, output_path, width, height, dpi, color_by_chrom):
    """Create piN vs piS scatter plot.

    \b
    Each gene is a point with piS on the x-axis and piN on the y-axis.
    Point size reflects gene length (codons).  A diagonal line marks
    piN = piS (neutral expectation).

    \b
    Examples:
      pie plot scatter -i gene_results.tsv -o scatter.png
      pie plot scatter -i gene_results.tsv -o scatter.png --color-by-chrom
    """
    from pie.plot import scatter_plot as _scatter_plot

    if not os.path.exists(input_path):
        click.echo(f"Error: file not found: {input_path}", err=True)
        sys.exit(1)

    _scatter_plot(input_path, output_path, width=width, height=height, dpi=dpi,
                  color_by_chrom=color_by_chrom)
    click.echo(f"Plot saved to {output_path}")


@plot.command(no_args_is_help=True, context_settings=_HELP_OPTS)
@_shared_plot_options
def histogram(input_path, output_path, width, height, dpi):
    """Create histogram of piN/piS distribution.

    \b
    Shows the distribution of piN/piS ratios across all genes, with a
    density curve overlay and a vertical line at piN/piS = 1.

    \b
    Examples:
      pie plot histogram -i gene_results.tsv -o histogram.png
    """
    from pie.plot import histogram_plot as _histogram_plot

    if not os.path.exists(input_path):
        click.echo(f"Error: file not found: {input_path}", err=True)
        sys.exit(1)

    _histogram_plot(input_path, output_path, width=width, height=height, dpi=dpi)
    click.echo(f"Plot saved to {output_path}")


@plot.command(no_args_is_help=True, context_settings=_HELP_OPTS)
@_shared_plot_options
def boxplot(input_path, output_path, width, height, dpi):
    """Create per-chromosome boxplots of piN, piS, and piN/piS.

    \b
    Three faceted panels showing the distribution of each metric per
    chromosome.  Useful for identifying chromosomes under different
    selection pressures.

    \b
    Examples:
      pie plot boxplot -i gene_results.tsv -o boxplot.png
    """
    from pie.plot import boxplot_plot as _boxplot_plot

    if not os.path.exists(input_path):
        click.echo(f"Error: file not found: {input_path}", err=True)
        sys.exit(1)

    _boxplot_plot(input_path, output_path, width=width, height=height, dpi=dpi)
    click.echo(f"Plot saved to {output_path}")


@plot.command("sliding-window", no_args_is_help=True, context_settings=_HELP_OPTS)
@_shared_plot_options
def sliding_window(input_path, output_path, width, height, dpi):
    """Create sliding window piN/piS line plot.

    \b
    Line plot of piN/piS along genomic coordinates, faceted by chromosome.
    Input should be window_results.tsv from 'pie run'.

    \b
    Examples:
      pie plot sliding-window -i window_results.tsv -o sw.png
    """
    from pie.plot import sliding_window_plot as _sliding_window_plot

    if not os.path.exists(input_path):
        click.echo(f"Error: file not found: {input_path}", err=True)
        sys.exit(1)

    _sliding_window_plot(input_path, output_path, width=width, height=height, dpi=dpi)
    click.echo(f"Plot saved to {output_path}")
```

**Step 4: Run tests to verify they pass**

Run: `pytest tests/test_cli.py::TestPlotCommand -v`
Expected: All PASS

**Step 5: Run all tests**

Run: `pytest tests/ -m "not slow" -v`
Expected: All PASS

**Step 6: Commit**

```bash
git add src/pie/cli.py tests/test_cli.py
git commit -m "Restructure pie plot as click.Group with 5 subcommands"
```

---

### Task 8: Final cleanup and full test run

**Files:**
- Verify: all modified files

**Step 1: Run all fast tests**

Run: `pytest tests/ -m "not slow" -v`
Expected: All PASS, no warnings about deprecated matplotlib usage

**Step 2: Verify CLI help renders correctly**

Run: `pie plot --help`
Expected: Shows all 5 subcommands

Run: `pie plot manhattan --help`
Expected: Shows shared + manhattan-specific options

**Step 3: Commit any remaining changes**

If any final adjustments were needed, commit them.
