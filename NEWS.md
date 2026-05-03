# cellchatr 0.1.0

## New Features

* `rust_subset_genes()` — HashSet-based gene filtering, O(n+m) vs R's O(n×m)
* `rust_wilcoxon_filter()` — parallel Wilcoxon rank-sum test with tie correction
  using Rayon, 183x faster than base R on 2,000 genes × 200 cells
* `rust_match_lr_pairs()` — HashSet-based ligand-receptor pair matching

## Validation

* All three functions validated on real `CellChatDB.human` (3,233 interactions)
* 16 unit tests — all passing
* R CMD check — zero errors, zero warnings
