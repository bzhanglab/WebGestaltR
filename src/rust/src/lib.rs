#![allow(non_snake_case)]
#![allow(clippy::too_many_arguments)]
use std::vec;

use ahash::{AHashMap, AHashSet};
use extendr_api::prelude::*;
use ndarray::Array2;
use webgestalt_lib::{
    methods::{
        gsea::{GSEAConfig, GSEAResult, RankListItem},
        multilist::{multilist_gsea, multilist_ora, GSEAJob, ORAJob},
        nta::{nta, NTAConfig},
        ora::{get_ora, ORAConfig, ORAResult},
    },
    readers::utils::Item,
};

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
/// @name fill_input_data_frame
/// @keywords internal
#[extendr]
pub fn fill_input_data_frame(gmt: Robj, genes: Robj, gene_sets: Robj) -> List {
    let genes_vec = genes.as_string_vector().unwrap();
    let gene_set_vec = gene_sets.as_string_vector().unwrap();
    let mut value_array = Array2::zeros((genes_vec.len(), gene_set_vec.len()));
    let mut gene_index: AHashMap<&String, usize> = AHashMap::default();
    let mut set_index: AHashMap<&String, usize> = AHashMap::default();
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
    for i in 0..value_array.len_of(ndarray::Axis(1)) {
        gene_set_val.push(value_array.column(i).to_vec())
    }
    // Construct DataFrame in R. Create list for now.
    List::from_names_and_values(gene_set_vec, gene_set_val).unwrap()
}

/// Calculate random walk permutations for a network from seeds
/// @name nta_rust
/// @param edge_list A list of edges
/// @param seeds A list of seeds
/// @return A list of nodes and scores
/// @author John Elizarraras
/// @keywords internal
#[extendr]
pub fn nta_rust(edge_list: Robj, seeds: Robj) -> List {
    println!("Running NTA Rust function.");
    let edge_list_as_vec: Vec<Vec<String>> = edge_list
        .as_list()
        .unwrap()
        .iter()
        .map(|(_, row)| row.as_string_vector().unwrap_or(vec![]))
        .filter(|x| !x.is_empty())
        .collect();
    println!("Creating Seeds");
    let seeds: Vec<String> = seeds.as_string_vector().unwrap();
    let config: NTAConfig = NTAConfig {
        edge_list: edge_list_as_vec,
        seeds,
        reset_probability: 0.5,
        tolerance: 1e-6,
        ..Default::default()
    };
    let res = nta(config);
    let nodes: Vec<String> = res.iter().map(|(x, _)| x.to_string()).collect();
    let scores: Vec<f64> = res.iter().map(|(_, x)| *x).collect();
    list!(nodes = nodes, scores = scores)
}

/// Run ORA using Rust library
/// @name ora_rust
/// @param sets A vector of analyte set names
/// @param parts A list of the analyte in the analyte sets
/// @param interest A vector of analytes of interest
/// @param reference A vector of analytes in the reference set
/// @returns A list of the results of ORA
/// @author John Elizarraras
/// @keywords internal
#[extendr]
fn ora_rust(sets: Robj, parts: Robj, interest: Robj, reference: Robj) -> List {
    let config: ORAConfig = ORAConfig {
        fdr_method: webgestalt_lib::stat::AdjustmentMethod::None,
        ..Default::default()
    };
    let mut gmt: Vec<Item> = Vec::new();
    let set_vec = sets.as_str_vector().unwrap();
    let parts_vec: Vec<Vec<String>> = parts
        .as_list()
        .unwrap()
        .iter()
        .map(|(_, x)| x.as_string_vector().unwrap())
        .collect();
    for (i, set) in set_vec.iter().enumerate() {
        gmt.push(Item {
            id: set.to_string(),
            url: String::default(),
            parts: parts_vec[i].clone(),
        })
    }
    let interest_set: AHashSet<String> = AHashSet::from_iter(interest.as_string_vector().unwrap());
    let reference_set: AHashSet<String> =
        AHashSet::from_iter(reference.as_string_vector().unwrap());
    let res: Vec<ORAResult> = get_ora(&interest_set, &reference_set, gmt, config);
    let mut p: Vec<f64> = Vec::new();
    let mut fdr: Vec<f64> = Vec::new();
    let mut expect: Vec<f64> = Vec::new();
    let mut enrichment_ratio: Vec<f64> = Vec::new();
    let mut overlap: Vec<i64> = Vec::new();
    let mut gene_set: Vec<String> = Vec::new();
    for row in res {
        gene_set.push(row.set);
        p.push(row.p);
        fdr.push(row.fdr);
        expect.push(row.expected);
        overlap.push(row.overlap);
        enrichment_ratio.push(row.enrichment_ratio);
    }
    list!(
        p = p,
        gene_set = gene_set,
        fdr = fdr,
        expect = expect,
        overlap = overlap,
        enrichment_ratio = enrichment_ratio
    )
}

/// Run multiomics ORA using Rust library
/// @param sets list of  the names of the analyte sets
/// @param big_part_vec list of the analyte in the analyte sets
/// @param interest list of analytes of interest
/// @param reference list of analytes in the reference set
/// @param method meta-analysis method to get meta-p values
/// @returns A list of vectors containing the results of ORA, with each list corresponding to each input list
/// @author John Elizarraras
/// @keywords internal
#[extendr]
pub fn rust_multiomics_ora(
    sets: Robj,
    big_part_vec: Robj,
    interest: Robj,
    reference: Robj,
    method: Robj,
) -> List {
    let config: ORAConfig = ORAConfig {
        fdr_method: webgestalt_lib::stat::AdjustmentMethod::None,
        ..Default::default()
    };
    let parts = big_part_vec.as_list().unwrap();
    let reference_lists = reference.as_list().unwrap();
    let method = match method.as_str().unwrap() {
        "fisher" => webgestalt_lib::methods::multilist::MultiListMethod::Meta(
            webgestalt_lib::methods::multilist::MetaAnalysisMethod::Fisher,
        ),
        _ => webgestalt_lib::methods::multilist::MultiListMethod::Meta(
            webgestalt_lib::methods::multilist::MetaAnalysisMethod::Stouffer,
        ),
    };
    let interest_vec = interest.as_list().unwrap();
    let big_set_vec = sets.as_list().unwrap();
    let mut jobs: Vec<ORAJob> = Vec::new();
    for (i, (_, big_set)) in big_set_vec.into_iter().enumerate() {
        let mut gmt: Vec<Item> = Vec::new();
        let set_vec = big_set.as_str_vector().unwrap();
        let parts_vec: Vec<Vec<String>> = parts[i]
            .as_list()
            .unwrap()
            .iter()
            .map(|(_, x)| x.as_string_vector().unwrap())
            .collect();
        for (i, set) in set_vec.iter().enumerate() {
            gmt.push(Item {
                id: set.to_string(),
                url: String::default(),
                parts: parts_vec[i].clone(),
            })
        }
        let interest_set: AHashSet<String> =
            AHashSet::from_iter(interest_vec[i].as_string_vector().unwrap());
        let reference_set: AHashSet<String> =
            AHashSet::from_iter(reference_lists[i].as_string_vector().unwrap());
        let job = ORAJob {
            gmt: gmt.clone(),
            interest_list: interest_set.clone(),
            reference_list: reference_set.clone(),
            config: config.clone(),
        };
        jobs.push(job)
    }
    let res: Vec<Vec<ORAResult>> =
        multilist_ora(jobs, method, webgestalt_lib::stat::AdjustmentMethod::None);
    let mut all_res: Vec<List> = Vec::new();
    for analysis in res {
        let mut p: Vec<f64> = Vec::new();
        let mut fdr: Vec<f64> = Vec::new();
        let mut expect: Vec<f64> = Vec::new();
        let mut enrichment_ratio: Vec<f64> = Vec::new();
        let mut overlap: Vec<i64> = Vec::new();
        let mut gene_set: Vec<String> = Vec::new();
        for row in analysis {
            gene_set.push(row.set);
            p.push(row.p);
            fdr.push(row.fdr);
            expect.push(row.expected);
            overlap.push(row.overlap);
            enrichment_ratio.push(row.enrichment_ratio);
        }
        all_res.push(list!(
            p = p,
            gene_set = gene_set,
            fdr = fdr,
            expect = expect,
            overlap = overlap,
            enrichment_ratio = enrichment_ratio
        ));
    }
    List::from_values(all_res)
}

/// Run GSEA using rust library
/// @param min_overlap the minimum overlap between analyte set and analyte list
/// @param max_overlap the maximum overlap between analyte set and analyte list
/// @param permutations the number of permutations to run
/// @param sets A vector of analyte set names
/// @param parts A list of the analytse in the analyte sets
/// @param analytes A vector of analytes names in the GSEA list
/// @param ranks A vector of ranks for the analytes in the GSEA list
/// @return List of the results of GSEA
/// @author John Elizarraras
/// @keywords internal
/// @name gsea_rust
#[extendr]
fn gsea_rust(
    min_overlap: Robj,
    max_overlap: Robj,
    permutations: Robj,
    sets: Robj,
    parts: Robj,
    analytes: Robj,
    ranks: Robj,
) -> List {
    let config = GSEAConfig {
        min_overlap: min_overlap.as_real().unwrap() as i32,
        max_overlap: max_overlap.as_real().unwrap() as i32,
        permutations: permutations.as_real().unwrap() as i32,
        ..Default::default()
    };
    let mut gmt: Vec<Item> = Vec::new();
    let set_vec = sets.as_str_vector().unwrap();
    let parts_vec: Vec<Vec<String>> = parts
        .as_list()
        .unwrap()
        .iter()
        .map(|(_, x)| x.as_string_vector().unwrap())
        .collect();
    for (i, set) in set_vec.iter().enumerate() {
        gmt.push(Item {
            id: set.to_string(),
            url: String::default(),
            parts: parts_vec[i].clone(),
        })
    }
    let mut analyte_list: Vec<RankListItem> = Vec::new();
    let analyte_vec: Vec<&str> = analytes.as_str_vector().unwrap();
    let ranks_vec: Vec<f64> = ranks.as_real_vector().unwrap();
    for (i, analyte) in analyte_vec.iter().enumerate() {
        analyte_list.push(RankListItem {
            analyte: analyte.to_string(),
            rank: ranks_vec[i],
        })
    }
    let res: Vec<GSEAResult> = webgestalt_lib::methods::gsea::gsea(analyte_list, gmt, config, None);
    let mut fdr: Vec<f64> = Vec::new();
    let mut p: Vec<f64> = Vec::new();
    let mut leading_edge: Vec<i32> = Vec::new();
    let mut gene_sets: Vec<String> = Vec::new();
    let mut es: Vec<f64> = Vec::new();
    let mut nes: Vec<f64> = Vec::new();
    let mut running_sum: Vec<Robj> = Vec::new();
    for row in res {
        fdr.push(row.fdr);
        p.push(row.p);
        leading_edge.push(row.leading_edge);
        gene_sets.push(row.set);
        es.push(row.es);
        nes.push(row.nes);
        running_sum.push(Robj::from(row.running_sum));
    }
    list!(
        fdr = fdr,
        p_val = p,
        ES = es,
        NES = nes,
        leading_edge = leading_edge,
        gene_sets = gene_sets.clone(),
        running_sum = List::from_names_and_values(gene_sets, running_sum),
    )
}

/// Run multiomics GSEA using rust library
/// @param min_overlap the minimum overlap between analyte set and analyte list
/// @param max_overlap the maximum overlap between analyte set and analyte list
/// @param permutations the number of permutations to run
/// @param sets A vector of analyte set names
/// @param parts A list of the analytse in the analyte sets
/// @param analytes A vector of analytes names in the GSEA list
/// @param ranks A vector of ranks for the analytes in the GSEA list
/// @param method_modifier method modifier for the multiomics method ("fisher" or "stouffer" if
/// meta-analysis, "mean", "max", or "rank" if other combination method)
/// @param combo_method method for combining analyte sets (meta, mean, or max)
/// @return List of lists of the results of GSEA
/// @author John Elizarraras
/// @keywords internal
/// @name gsea_rust
#[extendr]
pub fn rust_multiomics_gsea(
    min_overlap: Robj,
    max_overlap: Robj,
    permutations: Robj,
    sets: Robj,
    parts: Robj,
    analytes: Robj,
    ranks: Robj,
    method_modifier: Robj,
    combo_method: Robj,
) -> List {
    let config = GSEAConfig {
        min_overlap: min_overlap.as_real().unwrap() as i32,
        max_overlap: max_overlap.as_real().unwrap() as i32,
        permutations: permutations.as_real().unwrap() as i32,
        ..Default::default()
    };
    let method = if combo_method.as_str().unwrap() == "meta" {
        match method_modifier.as_str().unwrap() {
            "fisher" => webgestalt_lib::methods::multilist::MultiListMethod::Meta(
                webgestalt_lib::methods::multilist::MetaAnalysisMethod::Fisher,
            ),
            _ => webgestalt_lib::methods::multilist::MultiListMethod::Meta(
                webgestalt_lib::methods::multilist::MetaAnalysisMethod::Stouffer,
            ),
        }
    } else {
        let norm = match method_modifier.as_str().unwrap() {
            "mean" => webgestalt_lib::methods::multilist::NormalizationMethod::MeanValue,
            "median" => webgestalt_lib::methods::multilist::NormalizationMethod::MedianValue,
            "rank" => webgestalt_lib::methods::multilist::NormalizationMethod::MedianRank,
            _ => webgestalt_lib::methods::multilist::NormalizationMethod::None,
        };
        match combo_method.as_str().unwrap() {
            "max" => webgestalt_lib::methods::multilist::MultiListMethod::Max(norm),
            "mean" => webgestalt_lib::methods::multilist::MultiListMethod::Mean(norm),
            _ => webgestalt_lib::methods::multilist::MultiListMethod::Mean(norm),
        }
    };
    let mut gmt: Vec<Item> = Vec::new();
    let set_vec = sets.as_str_vector().unwrap();
    let parts_vec: Vec<Vec<String>> = parts
        .as_list()
        .unwrap()
        .iter()
        .map(|(_, x)| x.as_string_vector().unwrap())
        .collect();
    for (i, set) in set_vec.iter().enumerate() {
        gmt.push(Item {
            id: set.to_string(),
            url: String::default(),
            parts: parts_vec[i].clone(),
        });
    }
    let rank_list = ranks.as_list().unwrap();
    let mut jobs: Vec<GSEAJob> = Vec::new();
    for (i, (_, list)) in analytes.as_list().unwrap().into_iter().enumerate() {
        let mut analyte_list: Vec<RankListItem> = Vec::new();
        let analyte_vec: Vec<&str> = list.as_str_vector().unwrap();
        let ranks_vec: Vec<f64> = rank_list[i].as_real_vector().unwrap();
        for (i, analyte) in analyte_vec.iter().enumerate() {
            analyte_list.push(RankListItem {
                analyte: analyte.to_string(),
                rank: ranks_vec[i],
            })
        }
        let job = GSEAJob {
            gmt: gmt.clone(),
            rank_list: analyte_list,
            config: config.clone(),
        };
        jobs.push(job);
    }
    let res: Vec<Vec<GSEAResult>> =
        multilist_gsea(jobs, method, webgestalt_lib::stat::AdjustmentMethod::None);
    let mut all_res: Vec<List> = Vec::new();
    for analysis in res {
        let mut fdr: Vec<f64> = Vec::new();
        let mut p: Vec<f64> = Vec::new();
        let mut leading_edge: Vec<i32> = Vec::new();
        let mut gene_sets: Vec<String> = Vec::new();
        let mut es: Vec<f64> = Vec::new();
        let mut nes: Vec<f64> = Vec::new();
        let mut running_sum: Vec<Robj> = Vec::new();
        for row in analysis {
            fdr.push(row.fdr);
            p.push(row.p);
            leading_edge.push(row.leading_edge);
            gene_sets.push(row.set);
            es.push(row.es);
            nes.push(row.nes);
            running_sum.push(Robj::from(row.running_sum));
        }
        all_res.push(list!(
            fdr = fdr,
            p_val = p,
            ES = es,
            NES = nes,
            leading_edge = leading_edge,
            gene_sets = gene_sets.clone(),
            running_sum = List::from_names_and_values(gene_sets, running_sum),
        ))
    }
    List::from_values(all_res)
}

// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod WebGestaltR;
    fn fill_input_data_frame;
    fn gsea_rust;
    fn ora_rust;
    fn rust_multiomics_ora;
    fn nta_rust;
}
