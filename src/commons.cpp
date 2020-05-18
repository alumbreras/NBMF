#include <RcppArmadillo.h>
#include "commons.h"

using namespace Rcpp;
using namespace arma;

// Converts a matrix representation of Z (FxN)
// into a tensor representation (FxKxN)
// NA values in the Z matrix are converted to 
// indicator vectors will all positions set to 0.
// (no component selected)
void matrix_to_tensor(const arma::imat& Z, arma::icube& Z_tensor){
  int F = Z.n_rows;
  int N = Z.n_cols;
  int K = Z.max();
  Z_tensor.zeros();
  
  //arma::icube Z_tensor(F,K+1,N);
  Z_tensor.zeros();
  for(int n = 0; n<N; n++){
    for(int f = 0; f<F; f++){
      int k = Z(f,n);
      if(k != NA_INTEGER){
        Z_tensor(f,k,n) = 1;
      }
    }
  }
}

void vb_matrix_to_tensor(const arma::imat& Z, arma::cube& Z_tensor){
  int F = Z.n_rows;
  int N = Z.n_cols;
  int K = Z.max();
  
  //arma::icube Z_tensor(F,K+1,N);
  Z_tensor.zeros();
  for(int n = 0; n<N; n++){
    for(int f = 0; f<F; f++){
      int k = Z(f,n);
      if(k != NA_INTEGER){
        Z_tensor(f,k,n) = 1;
      }
    }
  }
}



// [[Rcpp::export]]
SEXP matrix_to_tensor_R(const arma::imat& Z, int Kmax){
  int F = Z.n_rows;
  int N = Z.n_cols;
  int K = Kmax;
  
  arma::icube Z_tensor(F,K,N);
  Z_tensor.zeros();
  for(int n = 0; n<N; n++){
    for(int f = 0; f<F; f++){
      int k = Z(f,n);
      if(k != NA_INTEGER){
        Z_tensor(f,k,n) = 1;
      }
    }
  }
  
  return wrap(Z_tensor);
}
