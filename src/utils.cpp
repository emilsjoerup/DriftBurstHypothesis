#include <RcppArmadillo.h>

//This code is taken from the repository: https://github.com/coatless/r-to-armadillo which is a collection of R functions written in C++

//[[Rcpp::export]]
arma::vec cfilter(arma::vec x, arma::vec filter, int sides = 2, bool circular = false)
{
  
  int nx = x.n_elem;
  int nf = filter.n_elem;
  int nshift;
  
  if(sides == NA_INTEGER || circular == NA_LOGICAL)  Rcpp::stop("invalid input");
  
  double z, tmp;
  
  if(sides == 2){
    nshift = nf /2;
  }
  else{
    nshift = 0;
  }
  
  arma::vec out = arma::zeros<arma::vec>(nx);
  
  if(!circular) {
    for(int i = 0; i < nx; i++) {
      z = 0;
      if(i + nshift - (nf - 1) < 0 || i + nshift >= nx) {
        out(i) = NA_REAL;
        continue;
      }
      for(int j = std::max(0, nshift + i - nx); j < std::min(nf, i + nshift + 1) ; j++) {
        tmp = x(i + nshift - j);
        z += filter(j) * tmp;
      }
      out(i) = z;
    }
  } else { /* circular */
    for(int i = 0; i < nx; i++)
    {
      z = 0;
      for(int j = 0; j < nf; j++) {
        int ii = i + nshift - j;
        if(ii < 0) ii += nx;
        if(ii >= nx) ii -= nx;
        tmp = x(ii);
        z += filter(j) * tmp;
      }
      out(i) = z;
    }
  }
  return out;
}
