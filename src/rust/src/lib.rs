use extendr_api::prelude::*;
use std::collections::HashSet;
use rayon::prelude::*;

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

/// Identify over-expressed genes per cell type using Wilcoxon rank-sum test.
/// @param counts NumericMatrix of gene expression (genes x cells).
/// @param labels Character vector of cell type labels, one per cell.
/// @param gene_names Character vector of gene names.
/// @param pval_threshold P-value cutoff for significance.
/// @return Character vector of over-expressed gene names.
/// @export
#[extendr]
fn rust_wilcoxon_filter(
    counts: RMatrix<f64>,
    labels: Vec<String>,
    gene_names: Vec<String>,
    pval_threshold: f64,
) -> Vec<String> {
    let ngenes = counts.nrows();
    let ncells = counts.ncols();

    // get unique cell types
    let cell_types: Vec<String> = {
        let mut seen = HashSet::new();
        labels.iter().filter(|l| seen.insert(*l)).cloned().collect()
    };

    // for each gene, test if it is overexpressed in any cell type
    let data: Vec<f64> = counts.data().to_vec();

    let overexpressed: Vec<String> = (0..ngenes)
        .into_par_iter()
        .filter_map(|gene_idx| {
            // extract expression values for this gene across all cells
            let gene_expr: Vec<f64> = (0..ncells)
                .map(|cell_idx| data[gene_idx + cell_idx * ngenes])
                .collect();

            // test against each cell type
            for ct in &cell_types {
                let group: Vec<f64> = labels.iter().zip(gene_expr.iter())
                    .filter(|(l, _)| *l == ct)
                    .map(|(_, v)| *v)
                    .collect();

                let other: Vec<f64> = labels.iter().zip(gene_expr.iter())
                    .filter(|(l, _)| *l != ct)
                    .map(|(_, v)| *v)
                    .collect();

                if group.is_empty() || other.is_empty() {
                    continue;
                }

                let pval = wilcoxon_pval(&group, &other);
                if pval < pval_threshold {
                    return Some(gene_names[gene_idx].clone());
                }
            }
            None
        })
        .collect();

    overexpressed
}

// Wilcoxon rank-sum test — returns approximate p-value using normal approximation
fn wilcoxon_pval(x: &[f64], y: &[f64]) -> f64 {
    let n1 = x.len() as f64;
    let n2 = y.len() as f64;

    // combine and rank
    let mut combined: Vec<(f64, usize)> = x.iter().map(|&v| (v, 0))
        .chain(y.iter().map(|&v| (v, 1)))
        .enumerate()
        .map(|(i, (v, g))| (v, g))
        .collect();

    combined.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

    // sum of ranks for group x (1-indexed)
    let w: f64 = combined.iter().enumerate()
        .filter(|(_, (_, g))| *g == 0)
        .map(|(rank, _)| rank as f64 + 1.0)
        .sum();

    let u = w - n1 * (n1 + 1.0) / 2.0;
    let mean_u = n1 * n2 / 2.0;
    let std_u = ((n1 * n2 * (n1 + n2 + 1.0)) / 12.0).sqrt();

    if std_u == 0.0 {
        return 1.0;
    }

    let z = (u - mean_u) / std_u;
    // two-tailed p-value using normal approximation
    2.0 * (1.0 - normal_cdf(z.abs()))
}

// Standard normal CDF approximation
fn normal_cdf(x: f64) -> f64 {
    0.5 * (1.0 + erf(x / std::f64::consts::SQRT_2))
}

fn erf(x: f64) -> f64 {
    let t = 1.0 / (1.0 + 0.3275911 * x.abs());
    let poly = t * (0.254829592
        + t * (-0.284496736
        + t * (1.421413741
        + t * (-1.453152027
        + t * 1.061405429))));
    let result = 1.0 - poly * (-x * x).exp();
    if x >= 0.0 { result } else { -result }
}

extendr_module! {
    mod cellchatr;
    fn hello_world;
    fn rust_subset_genes;
    fn rust_wilcoxon_filter;
}
