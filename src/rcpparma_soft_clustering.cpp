// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"


// [[Rcpp::depends(RcppArmadillo)]]

//' @title Use soft clustering algorithm
//
//' @param fr Weight vector
//' @param Zr The Z matrix
//' @param Sr the S matrix
//' @param nloops the number of loops to execute
//'
//' @return a list containing rhos and Zs
//' @export
//' @examples
//' n <- 3
//' power = -.5
//' M <- matrix(seq(9), nrow = n)
//' f <- rowSums(M) / sum(M)
//' D <- dist(f*M)
//' Delta <- as.numeric(0.5 * t(f) %*% D %*% f)
//' beta_rel <- 10^power
//' beta <- beta_rel * as.numeric(1 / Delta)
//' S <- as.matrix(exp(-beta * D))
//' index <- which.min(D %*% f)
//' Z <- matrix(0, n, n)
//' Z[, index] <- 1
//' Z <-  1e-20 * matrix(1, n, n) + (1 -  1e-20) * Z
//' soft_clustering(f, Z, S, Nloop=1000)
// [[Rcpp::export]]
Rcpp::List soft_clustering(Rcpp::NumericVector fr, Rcpp::NumericMatrix Zr, Rcpp::NumericMatrix Sr, int nloops) {

    int n = Sr.nrow();
    arma::mat S(Sr.begin(), n, n, false);
    arma::mat Z(Zr.begin(), n, n, false);
    arma::colvec f(fr.begin(), fr.size(), false);

    arma::colvec rho(n);
    arma::mat sxdiagr(n,n);

    for (int i = 0; i < nloops; i++) {
        rho = arma::trans(Z) * f;
        sxdiagr = S * arma::diagmat(rho);
        Z =  arma::diagmat(1 / arma::sum(sxdiagr,1)) * sxdiagr;
    }
    return Rcpp::List::create(
        Rcpp::Named("rho") = rho,
        Rcpp::Named("Z") = Z
    );
}
