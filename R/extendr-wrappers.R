#' Return string `"Hello world!"` to R.

# nolint start

#' @export
hello_world <- function() .Call(wrap__hello_world)

#' Filter expression matrix genes to only those present in CellChatDB.
#' @param expr_genes Character vector of gene names from expression matrix.
#' @param db_genes Character vector of gene names from CellChatDB.
#' @return Character vector of genes present in both inputs.
#' @export
rust_subset_genes <- function(expr_genes, db_genes) .Call(wrap__rust_subset_genes, expr_genes, db_genes)

#' Identify over-expressed genes per cell type using Wilcoxon rank-sum test.
#' @param counts NumericMatrix of gene expression (genes x cells).
#' @param labels Character vector of cell type labels, one per cell.
#' @param gene_names Character vector of gene names.
#' @param pval_threshold P-value cutoff for significance.
#' @return Character vector of over-expressed gene names.
#' @export
rust_wilcoxon_filter <- function(counts, labels, gene_names, pval_threshold) .Call(wrap__rust_wilcoxon_filter, counts, labels, gene_names, pval_threshold)


# nolint end
