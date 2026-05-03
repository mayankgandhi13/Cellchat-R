#' Return string `"Hello world!"` to R.

# nolint start

#' @return A character string "Hello world!"
#' @examples
#' hello_world()
#' @export
hello_world <- function() .Call(wrap__hello_world)

#' Filter expression matrix genes to only those present in CellChatDB.
#'
#' # Complexity
#' - Time:  O(n + m) where n = db_genes, m = expr_genes
#' - Space: O(n) for HashSet
#' - R equivalent (%%in%%): O(n x m)
#'
#' @param expr_genes Character vector of gene names from expression matrix.
#' @param db_genes Character vector of gene names from CellChatDB.
#' @return Character vector of genes present in both inputs.
#' @examples
#' rust_subset_genes(
#'   c("TGFB1", "VEGFA", "GAPDH"),
#'   c("TGFB1", "VEGFA", "MYC")
#' )
#' @export
rust_subset_genes <- function(expr_genes, db_genes) .Call(wrap__rust_subset_genes, expr_genes, db_genes)

#' Identify over-expressed genes per cell type using Wilcoxon rank-sum test.
#'
#' # Complexity
#' - Time:  O((g/p) x k x c log c) where g = genes, p = CPU cores,
#'          k = cell types, c = cells
#' - Space: O(c) per gene
#' - R equivalent: O(g x k x c log c) — serial, no parallelism
#'
#' @param counts NumericMatrix of gene expression (genes x cells).
#' @param labels Character vector of cell type labels, one per cell.
#' @param gene_names Character vector of gene names.
#' @param pval_threshold P-value cutoff for significance.
#' @return Character vector of over-expressed gene names.
#' @examples
#' counts <- matrix(c(rep(10, 50), rep(0, 50), rep(1, 100)),
#'   nrow = 2, byrow = TRUE)
#' labels <- c(rep("A", 50), rep("B", 50))
#' rust_wilcoxon_filter(counts, labels, c("TGFB1", "GAPDH"), 0.05)
#' @export
rust_wilcoxon_filter <- function(counts, labels, gene_names, pval_threshold) .Call(wrap__rust_wilcoxon_filter, counts, labels, gene_names, pval_threshold)

#' Match overexpressed genes against CellChatDB ligand-receptor pairs.
#'
#' # Complexity
#' - Time:  O(m + r) where m = overexpressed genes, r = LR pairs
#' - Space: O(m) for HashSet
#' - R equivalent (%%in%% twice): O(m x r)
#'
#' @param overexpressed Character vector of overexpressed gene names.
#' @param lr_ligands Character vector of ligands from CellChatDB.
#' @param lr_receptors Character vector of receptors from CellChatDB.
#' @param lr_names Character vector of interaction names from CellChatDB.
#' @return Character vector of matched interaction names.
#' @examples
#' rust_match_lr_pairs(
#'   c("TGFB1", "VEGFA"),
#'   c("TGFB1", "MYC", "VEGFA"),
#'   c("TGFBR1", "MYCBP", "FLT1"),
#'   c("TGFB1_TGFBR1", "MYC_MYCBP", "VEGFA_FLT1")
#' )
#' @export
rust_match_lr_pairs <- function(overexpressed, lr_ligands, lr_receptors, lr_names) .Call(wrap__rust_match_lr_pairs, overexpressed, lr_ligands, lr_receptors, lr_names)


# nolint end
