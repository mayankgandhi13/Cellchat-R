# cellchatr

An R package that accelerates [CellChat](https://github.com/jinworks/CellChat) 
database querying using Rust via [extendr](https://extendr.github.io/).

## Motivation

CellChat's gene filtering and statistical testing steps become serious 
bottlenecks on large single-cell datasets. This package reimplements three 
core functions in Rust — with parallel Wilcoxon testing via Rayon — and 
exposes them to R as drop-in replacements.

## Installation

```r
# install.packages("devtools")
devtools::install_github("mayankgandhi13/Cellchat-R", build_vignettes = TRUE)
```

## Quick Start

```r
library(cellchatr)
library(CellChat)

db <- CellChatDB.human

# Step 1: filter expression genes to CellChatDB
subset <- rust_subset_genes(expr_genes, db_genes)

# Step 2: identify overexpressed genes per cell type (parallel Wilcoxon)
overexpressed <- rust_wilcoxon_filter(
  counts         = your_matrix,
  labels         = cell_type_labels,
  gene_names     = subset,
  pval_threshold = 0.05
)

# Step 3: match to L-R pairs in CellChatDB
interactions <- rust_match_lr_pairs(
  overexpressed = overexpressed,
  lr_ligands    = db$interaction$ligand,
  lr_receptors  = db$interaction$receptor,
  lr_names      = db$interaction$interaction_name
)
```

See `vignette("cellchatr-workflow")` for a full worked example with 
real `CellChatDB.human` data.

## Functions

### `rust_subset_genes(expr_genes, db_genes)`
Filters expression matrix gene names to only those present in CellChatDB.

- **Algorithm:** Rust `HashSet` lookup
- **Complexity:** O(n + m) vs R's O(n × m)
- **Replaces:** `subsetData()`

### `rust_wilcoxon_filter(counts, labels, gene_names, pval_threshold)`
Runs a one-sided Wilcoxon rank-sum test with tie correction for each 
gene × cell type combination in parallel using Rayon.

- **Algorithm:** Parallel Wilcoxon with average-rank tie correction
- **Complexity:** O((g/p) × k × c log c) where p = CPU cores
- **Replaces:** `identifyOverExpressedGenes()`

### `rust_match_lr_pairs(overexpressed, lr_ligands, lr_receptors, lr_names)`
Matches overexpressed genes against CellChatDB ligand-receptor pairs.
An interaction is kept if its ligand OR receptor is overexpressed.

- **Algorithm:** Rust `HashSet` lookup
- **Complexity:** O(m + r) vs R's O(m × r)
- **Replaces:** `identifyOverExpressedInteractions()`

## Benchmarks

Tested on MacBook Air M-series (aarch64-apple-darwin, 8 cores):

### Wilcoxon filtering — the key bottleneck

| Dataset | Rust | Base R | Speedup |
|---|---|---|---|
| 2,000 genes × 200 cells × 2 cell types | 4.5ms | 822ms | **183x** |
| 2,000 genes × 200 cells × 2 cell types (vignette) | 57ms | 900ms | **16x** |
| 10,000 genes × 500 cells × 5 cell types | 100ms | ~5.5 hrs (est.) | **>1000x** |

### Gene lookup functions

| Function | Rust | Base R | Note |
|---|---|---|---|
| `rust_subset_genes()` | 1.9ms | 0.29ms | R faster — FFI overhead dominates |
| `rust_match_lr_pairs()` | 2.5ms | 0.38ms | R faster — FFI overhead dominates |

> **Why Rust wins on Wilcoxon:** 50,000 independent statistical tests 
> (10k genes × 5 cell types) run in parallel across all CPU cores. 
> R runs them serially one by one.

> **Why R wins on lookups:** R's internal C backend is already highly 
> optimized for simple operations. The Rust FFI crossing cost outweighs 
> the computation for small tasks.

## Validation

All three functions validated on real `CellChatDB.human` (3,233 L-R interactions):
- `rust_subset_genes()` — correctly filters 15,462 → 1,462 CellChatDB genes
- `rust_wilcoxon_filter()` — correctly detects TGFB1, TGFB2, TGFB3 as 
  overexpressed in simulated Fibroblasts
- `rust_match_lr_pairs()` — correctly maps 10 overexpressed genes to 89 
  interactions across 13 pathways (COLLAGEN, TGFb, VEGF, SPP1...)

## Tests

```r
devtools::test()
# [ FAIL 0 | WARN 0 | SKIP 0 | PASS 16 ]
```

16 edge case tests covering empty inputs, boundary thresholds, 
duplicate genes, and correct biological outputs.

## Requirements

- R >= 4.0
- Rust >= 1.65 (for building from source)
- CellChat >= 2.0 (for real database access)

## License

MIT
