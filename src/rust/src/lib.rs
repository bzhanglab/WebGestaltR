use extendr_api::prelude::*;
use extendr_api::wrapper::dataframe::Dataframe;
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
/// @export
#[extendr]
fn gsea_rust() -> () {
    // webgestalt_lib::methods::gsea::
}

/// Fill relation data frame for GSEA input
///
/// Fill 1 for gene in gene set
/// ## Diagram
/// ```shell
///         Gene Sets
///       ┌───────────┐  First column named 'gene' containing gene name
///       │A0100110100│  1 = in set
/// Genes │B0100101000│  0 = not in set
///       │C1011101001│
///       └───────────┘
/// ```
/// @param gmt A Data Frame with geneSet and gene columns from the GMT file
/// @param genes A vector of genes
/// @param geneSets A vector of gene sets
///
/// @return A Data Frame with the first column of gene and 1 or 0 for other columns of gene sets.
/// @author John Elizarraras
/// @keywords internal
/// @export
#[extendr]
pub fn fill_input_data_frame(gmt: Robj, genes: Robj, gene_sets: Robj) -> Robj {
    let genes_vec = genes.as_string_vector().unwrap();
    let mut gene_set_vec = gene_sets.as_string_vector().unwrap();
    let mut value_array = Array2::zeros((genes.len(), gmt.len()));
    let mut geneIndex: FxHashMap<&String, usize> = FxHashMap::default();
    let mut setIndex: FxHashMap<&String, usize> = FxHashMap::default();
    let gmtSet: Vec<String> = gmt.index("geneSet").unwrap().as_string_vector().unwrap();
    let gmtGene: Vec<String> = gmt.index("gene").unwrap().as_string_vector().unwrap();
    for (i, val) in genes_vec.iter().enumerate() {
        geneIndex.insert(val, i);
    }
    for (i, val) in gene_set_vec.iter().enumerate() {
        setIndex.insert(val, i);
    }
    for i in 0..gmtSet.len() {
        value_array[[geneIndex[&gmtGene[i]], setIndex[&gmtSet[i]]]] = 1;
    }
    gene_set_vec.insert(0, String::from("gene"));
    data_frame!(x = 1)
}

// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod WebGestaltR;
    fn rust_hello_world;
}
