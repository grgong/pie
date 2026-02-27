# Individual Sequencing Mode for pie

**Date:** 2026-02-27
**Status:** Approved

## Goal

Support multiple individual-sequencing samples by deriving pooled allele
frequencies from GT fields across selected samples, then feeding them into
the existing NG86 pipeline. This is mathematically equivalent to pairwise
population diversity estimation.

## CLI Changes

New `--mode` option on `pie run` (default `pool`):

```
pie run --mode pool ...         # existing behavior (default)
pie run --mode individual ...   # GT-based ("ind" accepted as alias)
```

### Mode-specific options

| Option           | Pool | Individual | Notes                              |
|------------------|------|------------|------------------------------------|
| `--sample`       | yes  | no         | Single sample for pool-seq         |
| `--min-depth`    | yes  | no         | Pool-seq read depth                |
| `--samples`      | no   | yes        | Comma-separated sample names       |
| `--samples-file` | no   | yes        | File, one sample per line          |
| `--min-call-rate`| no   | yes        | Default 0.8, range [0,1]           |
| `--min-an`       | no   | yes        | Minimum allele number, default 2   |

Shared options (both modes): `--vcf`, `--gff`, `--fasta`, `--outdir`,
`--min-freq`, `--min-qual`, `--pass-only`, `--keep-multiallelic`,
`--include-stop-codons`, `--window-size`, `--window-step`, `--threads`.

### Validation rules

1. `--samples` and `--samples-file` are mutually exclusive.
2. `--mode pool` + any individual-only option → error.
3. `--mode individual` + any pool-only option → error.
4. `--min-call-rate` must be in [0, 1].
5. `--mode individual` with neither `--samples` nor `--samples-file` → use
   all samples in VCF.
6. Named samples not found in VCF → error listing available samples.

## Architecture: Approach A — Separate `IndividualVariantReader`

New class in `vcf.py` alongside existing `VariantReader`. Both expose the
same `fetch(chrom, start, end) -> list[Variant]` interface. Downstream code
(`diversity.py`, `parallel.py`) is unchanged — it consumes `Variant` objects
regardless of origin.

Shared logic (QUAL/PASS filtering, multiallelic grouping) extracted into
helper functions to avoid duplication.

### Individual mode filter pipeline (per VCF record)

Single-pass GT extraction — no per-ALT rescanning:

1. Apply QUAL/PASS/REF filters.
2. **One pass** through selected samples' GT fields:
   - Count called samples → `call_rate = called / n_selected`.
   - Accumulate `ref_count` and per-ALT `alt_counts` from alleles.
   - `AN = sum of all counted alleles` (2 × called for diploid).
3. Check `call_rate >= --min-call-rate` → skip record if not.
4. Check `AN >= --min-an` → skip record if not.
5. For each ALT allele with `alt_count > 0`:
   - `freq = alt_count / AN`
   - Append `(pos, ref, alt, freq, AN, ref_count, alt_count)` to raw list.
6. Multiallelic grouping (same logic as pool mode).
7. Apply `--min-freq` filter at grouping phase.

`Variant.depth` is set to `AN` in individual mode (closest analogue to
pool-seq depth).

### Parallel runner changes

`run_parallel` gains: `mode`, `samples`, `min_call_rate`, `min_an`.
`_worker_init` instantiates `VariantReader` or `IndividualVariantReader`
based on mode into the same `_vcf` global. `_process_gene` and
`compute_gene_diversity` require no changes.

## Output Changes

### gene_results.tsv — individual mode adds:

| Column           | Type  | Description                                          |
|------------------|-------|------------------------------------------------------|
| `n_samples`      | int   | Number of selected samples                           |
| `mean_call_rate` | float | Mean call rate across variant sites in this gene     |

### summary.tsv — individual mode adds:

| Column               | Type  | Description                                          |
|----------------------|-------|------------------------------------------------------|
| `n_samples_selected` | int   | Total selected samples                               |
| `mean_call_rate`     | float | Variant-site-weighted mean call rate across all genes|

Pool mode: these columns are omitted to preserve backward compatibility.

### Implementation

`GeneResult` gets optional fields: `n_samples: int | None = None` and
`call_rates: list[float] | None = None` (per-variant-site call rates).
`mean_call_rate` is derived. Writers include the columns only when any
result has non-None values.

## Testing

1. **New fixtures** — multi-sample VCF with GT fields (3–4 diploid samples,
   some missing GTs).
2. **Unit tests for `IndividualVariantReader`:**
   - GT-to-frequency computation (hand-calculated)
   - `min_call_rate` filtering
   - `min_an` filtering
   - `min_freq` filtering on GT-derived AF
   - Missing genotype handling
   - `--samples` subset vs all samples
   - Multiallelic in individual mode
3. **CLI validation tests** — mode/option cross-validation errors.
4. **Integration test** — end-to-end `pie run --mode individual`.
5. **Pool mode regression** — all existing tests pass unchanged.
