#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//' EM algorithm
//'
//' Run EM using RcppArmadillo
//'
//' @param counts A diagonal matrix of TCC read counts
//' @param compatiblity_matrix A compatibility matrix between TCCs and transcripts
//' @param use_length_normalization A logical scalar specifying if normalization
//' using effective transcript lengths should be used during EM
//' @param tx_effective_lengths A numeric vector of effective transcript lengths
//' @param maxiter An integer scalar specifying the maximum number of EM iterations
//' @param conv A numeric scalar specifying the EM convergence threshold
//' @return A numeric vector of transcript read counts estimated using EM
//' @keywords internal
// [[Rcpp::export]]
arma::rowvec EM(const arma::mat counts,
                const arma::mat compatibility_matrix,
                const bool use_length_normalization,
                const arma::rowvec tx_effective_lengths,
                int maxiter,
                double conv) {

    arma::mat::elem_type counts_sum = accu(counts);
    arma::rowvec tx_fractions(compatibility_matrix.n_cols, arma::fill::ones);
    arma::rowvec tx_counts;
    arma::rowvec prev_tx_fractions(tx_fractions.n_elem, arma::fill::zeros);
    arma::rowvec prev_tx_counts(tx_fractions.n_elem, arma::fill::zeros);
    arma::mat probs;
    arma::vec probs_row_sums;
    arma::mat counts_distribution;
    int i;
    for (i = 0; i < maxiter; i++) {
        // Calculate probabilities by multiplying the columns of the match matrix
        // by isoform fractions, and then normalizing the matrix row-wise
        probs = compatibility_matrix;
        probs.each_row() %= tx_fractions;
        probs_row_sums = sum(probs, 1);
        probs.each_col() /= probs_row_sums;
        probs.replace(arma::datum::nan, 0);
        // Distribute reads across isoforms by probability
        counts_distribution = counts * probs;
        // Calculate isoform read counts & fractions
        tx_counts = sum(counts_distribution, 0);
        tx_fractions = tx_counts / counts_sum;
        if (use_length_normalization) {
            tx_fractions = tx_fractions / tx_effective_lengths;
            tx_fractions = tx_fractions / accu(tx_fractions);
        }
        // Stop if convergence is achieved
        if (approx_equal(tx_fractions, prev_tx_fractions, "absdiff", conv))
            break;
        prev_tx_fractions = tx_fractions;
        prev_tx_counts = tx_counts;
    }
    return tx_counts;
}
