use std::vec;

use extendr_api::prelude::*;
use ndarray::Array2;
use rustc_hash::FxHashMap;
use webgestalt_lib::methods::*;
/// Return string `"Hello world!"` to R.
/// @export
#[extendr]
fn rust_hello_world() -> &'static str {
    "Hello world!"
}

/// Run GSEA using rust library
/// @return List of the results of GSEA
/// @export
#[extendr]
fn gsea_rust() -> List {
    // webgestalt_lib::methods::gsea::
    list!(
        fdr = vec![0.01, 0.05, 0.1],
        leading_edge = vec![4, 6, 4],
        gene_sets = vec!["GO1", "GO2", "GO3"]
    )
}

/// Fill relation data frame for GSEA input
///
/// Fill 1 for gene in gene set
///
/// See https://github.com/extendr/extendr/issues/612 for how to export DataFrame
///
/// ## Diagram
/// ```shell
///         Gene Sets
///       ┌───────────┐  First column named 'gene' containing gene name
///       │A0100110100│  1 = in set
/// Genes │B0100101000│  0 = not in set
///       │C1011101001│  Due to limitiations with extendr-api v 0.6.0,
///       └───────────┘  function returns a list, and the R package will
///                      add the first 'gene' column
/// ```
/// @param gmt A Data Frame with geneSet and gene columns from the GMT file
/// @param genes A vector of genes
/// @param gene_sets A vector of gene sets
///
/// @return A Data Frame with the first column of gene and 1 or 0 for other columns of gene sets.
/// @author John Elizarraras
/// @keywords internal
/// @export
#[extendr]
pub fn fill_input_data_frame(gmt: Robj, genes: Robj, gene_sets: Robj) -> List {
    let genes_vec = genes.as_string_vector().unwrap();
    let gene_set_vec = gene_sets.as_string_vector().unwrap();
    let mut value_array = Array2::zeros((genes.len(), gmt.len()));
    let mut gene_index: FxHashMap<&String, usize> = FxHashMap::default();
    let mut set_index: FxHashMap<&String, usize> = FxHashMap::default();
    let gmt_set: Vec<String> = gmt.index("geneSet").unwrap().as_string_vector().unwrap();
    let gmt_gene: Vec<String> = gmt.index("gene").unwrap().as_string_vector().unwrap();
    for (i, val) in genes_vec.iter().enumerate() {
        gene_index.insert(val, i);
    }
    for (i, val) in gene_set_vec.iter().enumerate() {
        set_index.insert(val, i);
    }
    for i in 0..gmt_set.len() {
        value_array[[gene_index[&gmt_gene[i]], set_index[&gmt_set[i]]]] = 1;
    }
    let mut gene_set_val: Vec<Vec<i32>> = Vec::new();
    // gene_set_val.push(genes_vec.into_iter().map(|x| SafeTypes::String(x)).collect());
    for i in 0..gmt_set.len() {
        gene_set_val.push(value_array.column(i).to_vec())
    }
    // gene_set_vec.insert(0, String::from("gene"));
    // Construct DataFrame in R. Create list for now.
    List::from_names_and_values(gene_set_vec, gene_set_val).unwrap()
    // data_frame!(x = 1)
}

// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod WebGestaltR;
    fn rust_hello_world;
}
