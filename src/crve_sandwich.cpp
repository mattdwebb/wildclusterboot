// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>

// [[Rcpp::export]]
arma::mat crve_sandwich(arma::mat X, arma::mat bread, arma::vec uhat, arma::umat clusterby, arma::ivec comb, arma::uvec G) {

  int k = X.n_cols;
  int n = X.n_rows;
  int d = clusterby.n_cols;

  arma::mat inner_prod;
  arma::mat summation(k, k);

  arma::mat sandwich(k, k);
  arma::mat partial_sandwich(k, k);
  sandwich.zeros();
  int g;
  double u;

  for(int i=0; i < d; i++){

    inner_prod.set_size(k, G(i));
    inner_prod.zeros();

    for(int m=0; m < clusterby.n_rows; m++){

      g = clusterby(m,i);

      u = uhat(m);

      for(int k=0; k < X.n_cols; k++){

        inner_prod(k, g) += u*X(m, k);

      }
    }

    summation.zeros();

    for(int h=0; h < G(i); h++){

      for(int l=0; l < k; l++){

        for(int j=0; j < k; j++){

          summation(j,l) += inner_prod(j,h) * inner_prod(l,h);

        }
      }
    }

    partial_sandwich = bread * summation * bread * comb(i) * G(i)/(G(i)-1) * (n-1)/(n-k);
    sandwich += partial_sandwich;

  }

  if(d > 1){

    arma::cx_vec eigval;
    arma::cx_mat eigvec;

    arma::eig_gen(eigval, eigvec, sandwich);
    arma::vec realval = arma::conv_to<arma::vec>::from(eigval);
    arma::mat realvec = arma::conv_to<arma::mat>::from(eigvec);
    arma::mat eigdiag(realval.n_elem, realval.n_elem);
    eigdiag.zeros();

    //
    for(int i=0; i < realval.n_elem; i++){

      if(realval(i) < 0.0){

        eigdiag(i,i) = 0;

      } else {

        eigdiag(i,i) = realval(i);

      }

    }

    sandwich = realvec * eigdiag * realvec.t();

  }

  return sandwich;

}

// [[Rcpp::export]]
arma::mat eigen_fix_cpp(arma::mat sandwich) {

  arma::cx_vec eigval;
  arma::cx_mat eigvec;

  arma::eig_gen(eigval, eigvec, sandwich);

  arma::vec realval = arma::conv_to<arma::vec>::from(eigval);
  arma::mat realvec = arma::conv_to<arma::mat>::from(eigvec);
  arma::mat eigdiag(realval.n_elem, realval.n_elem);
  eigdiag.zeros();
  //
  for(int i=0; i < realval.n_elem; i++){

    if(realval(i) < 0.0){

      eigdiag(i,i) = 0;

    } else {

      eigdiag(i,i) = realval(i);

    }

  }

  sandwich = realvec * eigdiag * realvec.t();

  return realvec * eigdiag * realvec.t();

}
