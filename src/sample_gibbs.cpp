#include <RcppArmadillo.h>
//#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;
using namespace arma;
//TODO: si el W es caso todo 1, sera que el C es siempre 1 en todos los K?


// [[Rcpp::export]]
void printVecs(arma::irowvec v1, arma::irowvec v2, int K) {
  arma::imat m(2,K);
  m.row(0) = v1;
  m.row(1) = v2;
  m.print("Matrix m");
  
  // or equally well (where you could also print v1 and/or v2
  Rcpp::Rcout << "Matrix M\n" << m; 
}


// [[Rcpp::export]]
SEXP sample_gibbs_cpp(const arma::ivec& v_n, 
                        const arma::mat& W, 
                        arma::imat C,
                        double alpha = 1, int iter=100,
                        double burnin = 0.5){
  int K = W.n_cols;
  int F = W.n_rows;
  int nburnin = std::floor(iter*burnin);

  arma::mat    C_samples_mean(F, K); // Store E[C[f,k]]
  C_samples_mean.zeros();
  double j = 0;
  
  arma::irowvec     L(K);
  arma::rowvec      probs(K);
  arma::irowvec     ans(K);
  
  arma::mat likelihood_matrix(F,K);
  arma::rowvec prior(K);
  
  for(int f=0; f<F; f++){
    likelihood_matrix.row(f) = v_n[f] * W.row(f) + (1-v_n[f]) * (1-W.row(f));
    }

  L = arma::sum(C, 0); // tells how many counts in each k
  C_samples_mean = C_samples_mean + conv_to<mat>::from(C); //init sample
  j = j + 1;
  
  // Resample the topic assigment k of each count
  for(int i=1; i<iter; i++){
    for(int f=0; f<F; f++){
      
      // Skip empty entries in V
      //if(v_n[f] == 0){
      //  continue;
      //}
      
      // Remove current value
      C.row(f).zeros();
      L = arma::sum(C, 0);
    
      // Compute probability of choosing each component k 
      probs = (alpha + conv_to<rowvec>::from(L)) % likelihood_matrix.row(f);
      probs = probs/sum(probs);

      //Sample new k for the occurrence
      ans = ans.zeros();
      rmultinom(1, probs.begin(), K, ans.begin());
    }
    // After burning, store samples
    // (or sufficient statistics, if samples won't be needed)
    if(i >= nburnin){
      C_samples_mean = C_samples_mean + conv_to<mat>::from(C);; 
      j = j + 1.0;
    }
  }
  C_samples_mean = C_samples_mean/j;
  return Rcpp::wrap(C_samples_mean);
}
