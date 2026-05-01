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


# nolint end
