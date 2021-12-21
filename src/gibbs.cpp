#include <Rcpp.h>
using namespace Rcpp;

//' @title A Gibbs sampler using Rcpp
//' @description A Gibbs sampler using Rcpp
//' @param N the number of samples
//' @return a random sample of size \code{n}
//' @examples
//' \dontrun{
//' rnC <- gibbs_cpp(100)
//' rnC
//' }
//' @export
// [[Rcpp::export]]
NumericMatrix gibbs_cpp(int N) {
  int n = 100;
  int a = 30;
  int b = 60;
  NumericMatrix x( N, 2);
  NumericVector xt (2);
  
  for (int i = 1; i < N; i++) {
    xt[0] = x( i-1, 0 );
    xt[1] = x( i-1, 1 );
    xt[0] = rbinom(1, n, xt[1])[0];
    xt[1] = rbeta(1, xt[0] + a, n - xt[0] + b)[0];
    x( i, 0 ) = xt[0];
    x( i, 1 ) = xt[1];
  }
  
  return x;
}