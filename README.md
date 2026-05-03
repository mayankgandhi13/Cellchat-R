# cellchatr

An R package that accelerates [CellChat](https://github.com/sqjin/CellChat) 
database querying using Rust via [extendr](https://extendr.github.io/).

## Motivation

CellChat's gene filtering and statistical testing steps are bottlenecks on 
large single-cell datasets. This package reimplements three core functions 
in Rust, with parallelism via Rayon, and exposes them to R as drop-in 
replacements.

## Benchmark

### Small dataset (2,000 genes × 200 cells, 2 cell types)

| Function | Rust | Base R | Speedup |
|---|---|---|---|
| `rust_wilcoxon_filter()` | 4.5ms | 822ms | **183x faster** |
| `rust_subset_genes()` | 1.9ms | 0.29ms | R faster (FFI overhead) |
| `rust_match_lr_pairs()` | 2.5ms | 0.38ms | R faster (FFI overhead) |

### Large dataset (10,000 genes × 500 cells, 5 cell types)

| Function | Rust | Estimated R | Speedup |
|---|---|---|---|
| `rust_wilcoxon_filter()` | **100ms** | ~5.5 hours | **>1000x faster** |

> Rust wins decisively on computation-heavy parallel workloads (Wilcoxon testing).
> For simple lookups, R's internal C backend is faster than the Rust FFI crossing cost.
> Tested on MacBook Air M-series (aarch64-apple-darwin).

## Installation

```r
# install.packages("devtools")
devtools::install_github("mayankgandhi13/Cellchat-R")
```

## Usage

```r
library(cellchatr)

# 1. Filter expression matrix genes to CellChatDB genes
subset <- rust_subset_genes(expr_genes, db_genes)

# 2. Identify overexpressed genes per cell type (parallel Wilcoxon)
overexpressed <- rust_wilcoxon_filter(
  counts         = your_matrix,
  labels         = cell_type_labels,
  gene_names     = gene_names,
  pval_threshold = 0.05
)

# 3. Match overexpressed genes to L-R pairs in CellChatDB
interactions <- rust_match_lr_pairs(
  overexpressed = overexpressed,
  lr_ligands    = CellChatDB$interaction$ligand,
  lr_receptors  = CellChatDB$interaction$receptor,
  lr_names      = CellChatDB$interaction$interaction_name
)
```

## Functions

**`rust_subset_genes(expr_genes, db_genes)`**  
Filters expression matrix gene names to only those present in CellChatDB 
using a Rust HashSet for O(1) lookup.

**`rust_wilcoxon_filter(counts, labels, gene_names, pval_threshold)`**  
Runs a one-sided Wilcoxon rank-sum test for each gene × cell type combination 
in parallel using Rayon. Returns overexpressed gene names.

**`rust_match_lr_pairs(overexpressed, lr_ligands, lr_receptors, lr_names)`**  
Matches overexpressed genes against CellChatDB ligand-receptor pairs using 
a Rust HashSet. Returns interaction names where ligand OR receptor is 
overexpressed.

## Requirements

- R >= 4.0
- Rust >= 1.65 (for building from source)

## License

MIT
