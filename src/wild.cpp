// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>

// [[Rcpp::export]]
arma::mat bread_cpp(arma::mat X) {

  arma::mat bread = arma::inv(arma::trans(X) * X);

  return bread;

}

// [[Rcpp::export]]
arma::mat beta_cpp(arma::mat X, arma::mat bread, arma::mat y_wild) {

  arma::mat beta = bread * arma::trans(X) * y_wild;

  return beta;

}

// [[Rcpp::export]]
arma::mat y_weights(arma::vec fitted_data, arma::vec uhat, arma::mat boot_weights, arma::uvec bootind) {

  int b;
  arma::mat y_wild(uhat.n_elem, boot_weights.n_cols);

  for(int i=0; i < bootind.n_elem; i++){

    b = bootind(i);

    for(int j=0; j < boot_weights.n_cols; j++){

      y_wild(i,j) = boot_weights(b,j) * uhat(i) + fitted_data(i);

    }
  }

  return y_wild;

}
