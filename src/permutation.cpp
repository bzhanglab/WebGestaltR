#include <Rcpp.h>
#include <stdlib.h>
using namespace Rcpp;

const double SMALL_NUM = 0.000001;

NumericMatrix shuffleRows(NumericMatrix mat, IntegerVector index) {
	int nrow = mat.nrow(), ncol = mat.ncol();
	NumericMatrix m(nrow, ncol);
	for (int i = 0; i < nrow; i++) {
		m(i, _) = mat(index[i]-1, _);
	}
	return m;
}

void matrixMultiplyVectorByColumn(NumericMatrix mat, NumericVector vec) {
	for (int j = 0; j < mat.ncol(); j++) {
		mat(_, j) = mat(_, j) * vec;
	}
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
	int num_gene = inset_scores.nrow();
	int num_set = inset_scores.ncol();
	NumericMatrix rand_tot(num_gene, num_set);
	NumericVector rand_res(3 * num_set); //rand_min, rand_max, rand_best concatenated
	double rand_min, rand_max;
	
	IntegerVector rand_index = sample(num_gene, num_gene);
	NumericMatrix rand_inset_scores = shuffleRows(inset_scores, rand_index);
	NumericMatrix rand_outset_scores = shuffleRows(outset_scores, rand_index);
	matrixMultiplyVectorByColumn(rand_inset_scores, expression_value);
	NumericVector rand_set_tot = colSums(rand_inset_scores);
	
	for (int i = 0; i < num_gene; i++) {
		rand_inset_scores(i, _) = (rand_inset_scores(i, _) / (rand_set_tot + SMALL_NUM)) - rand_outset_scores(i, _);
	}

	for (int j = 0; j < num_set; j++) {
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
