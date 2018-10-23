#include <Rcpp.h>
#include <unordered_map>
using namespace Rcpp;

//' Fill relation data frame for GSEA input
//'
//' Fill 1 for gene in gene set
//'
//' @param gmt A Data Frame with geneSet and gene columns from the GMT file
//' @param genes A vector of genes
//' @param geneSets A vector of gene sets
//'
//' @return A Data Frame with the first column of gene and 1 or 0 for other columns of gene sets.
//' @author Yuxing Liao
//' @keywords internal
// [[Rcpp::export]]
DataFrame fillInputDataFrame(DataFrame gmt, CharacterVector genes, CharacterVector geneSets) {
 	IntegerMatrix rel(genes.size(), geneSets.size());
	std::unordered_map<String, size_t> geneIndex, setIndex;
	CharacterVector gmtSet = as<CharacterVector>(gmt["geneSet"]);
	CharacterVector gmtGene = as<CharacterVector>(gmt["gene"]);
	List result(geneSets.size() + 1);

	for (size_t i = 0, length = genes.size(); i < length; i++) {
		geneIndex[genes[i]] = i;
	}
	for (size_t i = 0, length = geneSets.size(); i < length; i++) {
		setIndex[geneSets[i]] = i;
	}
	for (size_t i = 0, length = gmtSet.size(); i < length; i++) {
		rel(geneIndex[gmtGene[i]], setIndex[gmtSet[i]]) = 1;
	}
	result[0] = genes;
	for (size_t i = 1, length = result.size(); i < length; i++) {
		result[i] = rel(_, i - 1);
	}
	geneSets.push_front("gene"); //modify in place and back in R
	// Set colnames after creating DF. Use list names cannot avoid converting colon to
	// syntactic names (no 'optional' parameter like in as.data.frame for list)
	DataFrame df = DataFrame::create(result, _["stringsAsFactors"]=false);
	df.attr("names")= geneSets;
	return df;
}
