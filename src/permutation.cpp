#include <Rcpp.h>
#include <stdlib.h>
using namespace Rcpp;

const double SMALL_NUM = 0.000001;

NumericMatrix shuffleAndMultiplyColumn(NumericMatrix mat, NumericVector vec, IntegerVector rand_index) {
	//multiply vector with matrix columns which are shuffled ad hoc by index
	size_t nrow = mat.nrow(), ncol = mat.ncol();
	NumericMatrix m(nrow, ncol);
	for (size_t j = 0; j < ncol; j++) {
		for (size_t i = 0; i < nrow; i++) {
			m(i, j) = mat(rand_index[i] - 1, j) * vec[i];
		}
	}
	return m;
}

//' Permutaion in GSEA algorithm
//'
//' @param inset_scores Scaled score matrix for genes in sets
//' @param outset_scores Normalized score matrix for genes not in sets
//' @param expression_value Vector of gene rank scores
//'
//' @return A vector of concatenated random minimal,maimum and best running sum scores for each set.
//' @author Yuxing Liao
//' @keywords internal
// [[Rcpp::export]]
NumericVector gseaPermutation(NumericMatrix inset_scores, NumericMatrix outset_scores, NumericVector expression_value) {
	size_t num_gene = inset_scores.nrow();
	size_t num_set = inset_scores.ncol();
	NumericMatrix rand_tot(num_gene, num_set);
	NumericVector rand_res(3 * num_set); //rand_min, rand_max, rand_best concatenated
	double rand_min, rand_max;

	IntegerVector rand_index = sample(num_gene, num_gene);
	NumericMatrix rand_inset_scores = shuffleAndMultiplyColumn(inset_scores, expression_value, rand_index);
	NumericVector rand_set_tot = colSums(rand_inset_scores);

	for (size_t i = 0; i < num_gene; i++) {
		//shuffle outset scores here when in use without copying values beforehand
		rand_inset_scores(i, _) = (rand_inset_scores(i, _) / (rand_set_tot + SMALL_NUM)) - outset_scores(rand_index[i] - 1, _);
	}

	for (size_t j = 0; j < num_set; j++) {
		//explicit type conversion needed https://thecoatlessprofessor.com/programming/unofficial-rcpp-api-documentation/#carth
		rand_tot(_, j) = cumsum(rand_inset_scores(_, j)).get();
		rand_max = max(rand_tot(_, j));
		rand_min = min(rand_tot(_, j));
		if (rand_max > abs(rand_min)) {
			rand_res[j + 2 * num_set] = rand_max;
		} else {
			rand_res[j + 2 * num_set] = rand_min;
		}
		rand_res[j] = rand_min;
		rand_res[j + num_set] = rand_max;
	}
	return rand_res;
}
