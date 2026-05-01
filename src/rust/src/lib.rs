use extendr_api::prelude::*;
use std::collections::HashSet;

/// Return string `"Hello world!"` to R.
/// @export
#[extendr]
fn hello_world() -> &'static str {
    "Hello world!"
}

/// Filter expression matrix genes to only those present in CellChatDB.
/// @param expr_genes Character vector of gene names from expression matrix.
/// @param db_genes Character vector of gene names from CellChatDB.
/// @return Character vector of genes present in both inputs.
/// @export
#[extendr]
fn rust_subset_genes(expr_genes: Vec<String>, db_genes: Vec<String>) -> Vec<String> {
    let db_set: HashSet<String> = db_genes.into_iter().collect();
    expr_genes
        .into_iter()
        .filter(|g| db_set.contains(g))
        .collect()
}

extendr_module! {
    mod cellchatr;
    fn hello_world;
    fn rust_subset_genes;
}
