#include <RcppArmadillo.h>
#include "commons.h"
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/gamma.hpp>

//using digamma = boost::ios_base::fmtflags;
//using digamma = boost::math::special_functions::digamma;
using namespace Rcpp;
using namespace arma;


struct vbcounts {
   arma::mat f_ones; //declared as double to be ready for 'double' probabilities
   arma::mat f_zeros;
   arma::mat f_total;
   arma::mat n_total;
   arma::rowvec fn_total;
   int Kmax;
   
   explicit vbcounts(int f, int k, int n):  
                                         f_ones(n,k),
                                         f_zeros(n,k),
                                         f_total(n,k),
                                         n_total(f,k),
                                         fn_total(k),
                                         Kmax(k) {
                                                   f_ones.zeros();
                                                   f_zeros.zeros();
                                                   f_total.zeros();
                                                   n_total.zeros();
                                                   fn_total.zeros();
                                                 }
};


// Counts when Z is a FxN matrix
// Z: a FxN matrix that contains z_fn indicator vectors of size K
// Ignore NA's. NA are considered hidden or unknown elements,
// for instance to be used in testing predictions.
void compute_vbcounts(const arma::imat& V, 
                      const arma::cube& E_Z, 
                      vbcounts& cnt){

  cnt.f_ones.zeros();
  cnt.f_zeros.zeros();
  cnt.f_total.zeros();
  cnt.n_total.zeros();
  cnt.fn_total.zeros();
  
  int N = V.n_cols;
  int F = V.n_rows;
  for(int n=0; n<N; n++){
    for(int f=0; f<F; f++){
      if(V(f,n) != NA_INTEGER){
        arma::rowvec probs      = E_Z.slice(n).row(f);
        cnt.f_ones.row(n)      += probs * (V(f,n));
        cnt.f_zeros.row(n)     += probs * (1-V(f,n));
        cnt.f_total.row(n)     += probs;
        cnt.n_total.row(f)     += probs;
        cnt.fn_total           += probs;
      }
    }
  }
}

// Assign k to a (f,n) pair and update counts
void set_vbassignment(int f, int n, arma::rowvec probs, 
                    arma::cube& E_Z,
                    const arma::imat& V, vbcounts& cnt){
  E_Z.slice(n).row(f)   = probs;     
  cnt.f_ones.row(n)    += probs * V(f,n);
  cnt.f_zeros.row(n)   += probs * (1-V(f,n));
  cnt.f_total.row(n)   += probs;
  cnt.n_total.row(f)   += probs;
  cnt.fn_total         += probs;
}

// Remove assigment from (f,n) pair and update counts
void remove_vbassignment(int f, int n, 
                       arma::cube& E_Z,
                       const arma::imat& V, vbcounts& cnt){

  arma::rowvec probs       = E_Z.slice(n).row(f);  // current assignment
  E_Z.slice(n).row(f).zeros();
  cnt.f_ones.row(n)      -= probs * V(f,n);
  cnt.f_zeros.row(n)     -= probs * (1-V(f,n));
  cnt.f_total.row(n)     -= probs;
  cnt.n_total.row(f)     -= probs;
  cnt.fn_total           -= probs;
}



// [[Rcpp::export]]
SEXP VB_DirBer_Rcpp(const arma::imat& V, 
                              arma::imat Z, int K,
                              double alpha = 1, double beta =1, double gamma=1, 
                              int iter=100){
  int Kmax = K+1;
  int F = V.n_rows;
  int N = V.n_cols;
  arma::rowvec	probs(Kmax);
  arma::mat	    gamma_vb(F,Kmax);
  arma::mat 	  alpha_vb(N,Kmax);
  arma::mat   	beta_vb(N,Kmax);

  //Initial E_Z tensor from Z assignments matrix
  arma::cube E_Z(F,Kmax,N);
  vb_matrix_to_tensor(Z, E_Z);

  // First counts
  // While sampling we only update individual elements
  vbcounts cnt(F, Kmax, N);
  compute_vbcounts(V, E_Z, cnt);

  // Resample the topic assigment k of each count
  for(int i=1; i<iter; i++){
    Rcout << "iter: " << i << '/' << iter << std::endl;
    
    for(int n=0; n<N; n++){
      for(int f=0; f<F; f++){
        if(V(f,n) != NA_INTEGER){
          remove_vbassignment(f, n, E_Z, V, cnt);
          
          // Compute probability of every component
          gamma_vb.row(f) = gamma  + cnt.n_total.row(f);
          alpha_vb.row(n) = alpha  + cnt.f_ones.row(n);
          beta_vb.row(n)  = beta   + cnt.f_zeros.row(n);
          probs = gamma_vb.row(f) % 
                     (V(f,n) * alpha_vb.row(n) + (1-V(f,n)) * beta_vb.row(n)) /
                     (alpha_vb.row(n) + beta_vb.row(n));
          probs = probs/accu(probs);
          
          // Update E_Z(f,,n)
          set_vbassignment(f, n, probs, E_Z, V, cnt);
        }
      }
    }
  }
  
  // Assertions
  if(arma::any(arma::vectorise(alpha_vb) < 0)) {stop("Negative values in alpha_vb");}
  if(arma::any(arma::vectorise(beta_vb) < 0)) {stop("Negative values in beta_vb");}
  
  
  // Compute and return the expected factors
  // E_W: Expectation of a Multinomial
  // E_H: Expectation of a Beta
  arma::mat E_W(F,Kmax);
  arma::mat E_H(N,Kmax);
  
  for(int f=0; f<F; f++){
    E_W.row(f) = gamma_vb.row(f) / accu(gamma_vb.row(f)); 
  }
  for(int n=0; n<N; n++){
    E_H.row(n) = alpha_vb.row(n)  / (alpha_vb.row(n) + beta_vb.row(n));
  }
  
  // Assertions
  if(arma::any(arma::vectorise(E_W) < 0)) {stop("Negative values in E_W");}
  if(arma::any(arma::vectorise(E_H) < 0)) {stop("Negative values in E_H");}
  
  return Rcpp::List::create(Rcpp::Named("E_Z") = wrap(E_Z),
                            Rcpp::Named("E_W") = wrap(E_W),
                            Rcpp::Named("E_H") = wrap(E_H));
}



// [[Rcpp::export]]
double E_log_p_W_Rcpp(const arma::mat E_log_W, double gamma){
  int F = E_log_W.n_rows;
  int K = E_log_W.n_cols;
  double x = 0;
  x = F*(boost::math::lgamma(K*gamma) - K*boost::math::lgamma(gamma));
  for(int f=0; f < F; f++){
    x += (gamma-1) * arma::accu(E_log_W.row(f));
  }
  return(x);
}

// [[Rcpp::export]]
double E_log_q_W_Rcpp(const arma::mat& E_log_W, const arma::mat& gamma_vb){
  int F = E_log_W.n_rows;
  int K = E_log_W.n_cols;
  double x = 0;
  for(int f=0; f<F; f++){
    x += boost::math::lgamma(arma::accu(gamma_vb.row(f)));
    for(int k=0; k<K; k++){
      x -= boost::math::lgamma(gamma_vb(f,k));
    }
    x += arma::accu((gamma_vb.row(f)-1) % E_log_W.row(f));
  }
  return(x);
} 


// [[Rcpp::export]]
double E_log_p_Z_Rcpp(const arma::cube& E_Z, const arma::mat& E_log_W,
                      const arma::imat& V){
  int N = E_Z.n_slices;
  int F = E_Z.n_rows;
  int K = E_Z.n_cols;
  double x = 0;
  for(int n=0; n<N; n++){
    for(int f=0; f<F; f++){
      if(V(f,n) != NA_INTEGER){
        for(int k=0; k<K; k++){
          x += E_Z(f,k,n) * E_log_W(f,k);
          //Rcout << "(f,n):" << f << "," << n << ":" << E_Z(f,k,n) << std::endl;
          //Rcout << "E_Z(f,k,n):" << E_Z(f,k,n) << std::endl;
          //Rcout << "log E_Z(f,k,n):" << E_log_W(f,k) << std::endl;
          //Rcout << "x:" << x << std::endl;
          
        }
      }
    }
  }
  return(x);
}

// [[Rcpp::export]]
double E_log_q_Z_Rcpp(const arma::cube& E_Z, const arma::imat& V){
  int N = E_Z.n_slices;
  int F = E_Z.n_rows;
  int K = E_Z.n_cols;
  double x = 0;
  for(int n=0; n<N; n++){
    for(int f=0; f<F; f++){
      if(V(f,n) != NA_INTEGER){
        for(int k=0; k<K; k++){
          x += E_Z(f,k,n) * log(E_Z(f,k,n));
        }
      }
    }
  }
  return(x);
}

// [[Rcpp::export]]
double E_log_p_H_Rcpp(const arma::mat& E_log_H, const arma::mat& E_log_1_H, 
                     double alpha, double beta){
  
  int N = E_log_H.n_rows;
  int K = E_log_H.n_cols;
  double x = 0;
  x = N*K*(boost::math::lgamma(alpha+beta) - 
    (boost::math::lgamma(alpha) + boost::math::lgamma(beta))) +
    (alpha-1)*arma::accu(E_log_H) + (beta-1)*arma::accu(E_log_1_H);
  return(x);
}

// [[Rcpp::export]]                       
double E_log_q_H_Rcpp(const arma::mat& E_log_H, const arma::mat& E_log_1_H, 
                     const arma::mat& alpha_vb, const arma::mat& beta_vb){
  
  int N = E_log_H.n_rows;
  int K = E_log_H.n_cols;
  double x = 0;
  for(int n=0; n<N; n++){
    for(int k=0; k<K; k++){
      x += boost::math::lgamma(alpha_vb(n,k) + beta_vb(n,k)) -
          (boost::math::lgamma(alpha_vb(n,k)) + 
           boost::math::lgamma(beta_vb(n,k))) +
           (alpha_vb(n,k)-1) * E_log_H(n,k) + 
           (beta_vb(n,k)-1) * E_log_1_H(n,k);
    }
  }
  //x += arma::accu((alpha_vb-1) % E_log_H + (beta_vb-1) % E_log_1_H);
  return (x);
}

// [[Rcpp::export]]     
double E_log_p_V_Rcpp(const arma::imat& V, const arma::cube& E_Z, 
                     const arma::mat& E_log_H, const arma::mat& E_log_1_H){
  int N = V.n_cols;
  int F = V.n_rows;
  int K = E_log_H.n_cols;
  
  double x = 0;
  for(int n=0; n<N; n++){
    for(int f=0; f<F; f++){
      for(int k=0; k<K; k++){
        if(V(f,n) != NA_INTEGER){
          x += E_Z(f,k,n) * (V(f,n) * E_log_H(n,k) +
                            (1-V(f,n)) * E_log_1_H(n,k));
        }
      }
    }
  }
  return(x);
}

// [[Rcpp::export]]
double likelihood_aspect(const arma::imat& V,
                         const arma::mat& gamma_vb, 
                         const arma::mat& alpha_vb, 
                         const arma::mat& beta_vb){
  
  int F = V.n_rows;
  int N = V.n_cols;
  int K = gamma_vb.n_cols;
  
  double likelihood = 0;
  arma::mat E_WH;
  
  // Compute expectations
  // E_W: Expectation of a Multinomial
  // E_H: Expectation of a Beta
  arma::mat E_W(F,K);
  arma::mat E_H(N,K);
  for(int f=0; f<F; f++){
    E_W.row(f) = gamma_vb.row(f) / arma::accu(gamma_vb.row(f)); 
  }
  for(int n=0; n<N; n++){
    E_H.row(n) = alpha_vb.row(n)  / (alpha_vb.row(n) + beta_vb.row(n));
  }
  
  // Compute p(V) = \prod_{fn} v_{fn}*(E[WH]) + (1-v_{fn})(1-E[WH])
  E_WH = E_W * E_H.t();
  
  for(int f=0; f<F; f++){
    for(int n=0; n<N; n++){
      if(V(f,n) != NA_INTEGER){
        for(int k=0; k<K; k++){
          likelihood += V(f,n) * log(E_WH(f,n)) + 
                        (1-V(f,n)) * log(1-E_WH(f,n));
        }  
      }
    }
  }
  //likelihood = arma::accu(V % log(E_WH) + (1-V) % log(1-E_WH));
  return(likelihood);
}

// [[Rcpp::export]]
double lower_bound_aspect(const arma::cube& E_Z, 
                         const arma::mat& E_log_W,
                         const arma::mat& E_log_H,
                         const arma::mat& E_log_1_H,
                         const arma::mat& gamma_vb, 
                         const arma::mat& alpha_vb, const arma::mat& beta_vb, 
                         double alpha, double beta, double gamma, 
                         const arma::imat& V){
  //int Kmax = K+1;
  int F = V.n_rows;
  int N = V.n_cols;
  int K = E_log_W.n_cols;
  double E_log_p_W = 0;
  double E_log_p_Z = 0;
  double E_log_p_H = 0;
  double E_log_p_V = 0;
  double E_log_q_W = 0;
  double E_log_q_H = 0;
  double E_log_q_Z = 0;
  double lowerbound = 0;
  //Rcout << "computing lower bound:" << std::endl;
  
  E_log_p_W = E_log_p_W_Rcpp(E_log_W, gamma);
  E_log_q_W = E_log_q_W_Rcpp(E_log_W, gamma_vb);
  E_log_p_Z = E_log_p_Z_Rcpp(E_Z, E_log_W, V);
  E_log_q_Z = E_log_q_Z_Rcpp(E_Z, V);
  E_log_p_H = E_log_p_H_Rcpp(E_log_H, E_log_1_H, alpha, beta);
  E_log_q_H = E_log_q_H_Rcpp(E_log_H, E_log_1_H, alpha_vb, beta_vb);
  E_log_p_V = E_log_p_V_Rcpp(V, E_Z, E_log_H, E_log_1_H);

  lowerbound = E_log_p_W + E_log_p_Z + E_log_p_H + E_log_p_V -
               E_log_q_W - E_log_q_Z - E_log_q_H;
  
  return(lowerbound);
}
  

// [[Rcpp::export]]
SEXP VB_Aspect_Rcpp(const arma::imat& V, 
                    arma::imat Z, int K,
                    double alpha = 1, double beta =1, double gamma=1, 
                    int iter=100, bool checklb = false){
  
  //int Kmax = K+1;
  int Kmax = K;
  int F = V.n_rows;
  int N = V.n_cols;
  
  // Shouldn't lowerbounds be the same as likelihood p(V)?
  // Likelihood is computed as:
  // p(V) = v_{fn}*E[WH] + (1_v{fn})*(1-E[WH])
  arma::rowvec  lowerbounds(iter); // lower bound traces
  arma::rowvec  likelihoods(iter); // likelihood traces
  arma::rowvec  kldivergences(iter); // likelihood traces
  lowerbounds.fill(datum::nan);
  likelihoods.fill(datum::nan);
  kldivergences.fill(datum::nan);
  
  arma::rowvec	probs(Kmax);
  arma::mat	    gamma_vb(F,Kmax);
  arma::mat 	  alpha_vb(N,Kmax);
  arma::mat   	beta_vb(N,Kmax);
  arma::mat     E_log_W(F, Kmax);
  arma::mat     E_log_H(N, Kmax);
  arma::mat     E_log_1_H(N, Kmax);
  double lb = 0;
  double like = 0;
  float last_lb = -1000000;
  float last_like = -1000000;
    
  //Initial E_Z tensor from Z assignments matrix
  arma::cube E_Z(F,Kmax,N);
  vb_matrix_to_tensor(Z, E_Z);

  // First counts
  // While sampling we only update individual elements
  vbcounts cnt(F, Kmax, N);
  compute_vbcounts(V, E_Z, cnt);

  // Resample the topic assigment k of each count
  for(int i=0; i<iter; i++){
    Rcout << "iter: " << i << '/' << iter << std::endl;
    
    
    for(int f=0; f<F; f++){
      gamma_vb.row(f) = gamma  + cnt.n_total.row(f);
    }
    for(int n=0; n<N; n++){
      alpha_vb.row(n) = alpha  + cnt.f_ones.row(n);
      beta_vb.row(n)  = beta   + cnt.f_zeros.row(n);
    }
    
    // q(W)
    for(int f=0; f<F; f++){
      for(int k=0; k<Kmax; k++){
        E_log_W(f,k) =  boost::math::digamma(gamma_vb(f,k));
      }
      E_log_W.row(f) = E_log_W.row(f) - 
                        boost::math::digamma(arma::accu(gamma_vb.row(f)));
    }

    // q(H)
    for(int n=0; n<N; n++){
      for(int k=0; k<Kmax; k++){
        E_log_H(n,k)   =  boost::math::digamma(alpha_vb(n,k)) - 
                          boost::math::digamma(alpha_vb(n,k) + beta_vb(n,k));
        E_log_1_H(n,k) =  boost::math::digamma(beta_vb(n,k)) -
                          boost::math::digamma(alpha_vb(n,k) + beta_vb(n,k));
      }
    }
      
    // q(Z)
    for(int n=0; n<N; n++){
      for(int f=0; f<F; f++){
        if(V(f,n) != NA_INTEGER){
          remove_vbassignment(f, n, E_Z, V, cnt);

          // Compute expectation of every component
          probs = E_log_W.row(f) + 
                                    V(f,n)  * E_log_H.row(n) + 
                                    (1-V(f,n)) * E_log_1_H.row(n);
          probs = arma::exp(probs);
          probs = probs / arma::accu(probs);

          // Update E_Z(f,,n) and the counters
          set_vbassignment(f, n, probs, E_Z, V, cnt);
        } else{
          E_Z.slice(n).row(f).zeros(); // works as NA, to contribute 0 in the sums
        }
      }
    }
    /////////////////////////////////////////////////////////////
    // Monitorize lower bound (and likelihood)
    /////////////////////////////////////////////////////////////
    if(checklb){
      lb   = lower_bound_aspect(E_Z, E_log_W, E_log_H, E_log_1_H,
                                gamma_vb, alpha_vb, beta_vb, 
                                alpha, beta, gamma, V);
      Rcout << "lower bound: " << lb   << std::endl;
    } else {
      lb = 0;
    }
    lowerbounds[i] = lb; // lower bound traces
    
    
    // Likelihood behaves weirdly
    //like = likelihood_aspect(V, gamma_vb, alpha_vb, beta_vb);
    //likelihoods[i] = like;
    //kldivergences[i] = like - lb;
    //Rcout << "likelihood: "  << like << std::endl;
    //Rcout << "KL divergence: "  << like - lb << std::endl;
    
    if(i>1){
      if(lb < last_lb){
        Rcout << 
            "Decrease detected in lower bound. Something is wrong. Iter: " 
            << i << ". Decrease: " << lb - last_lb << std::endl;
        warning("Decrease detected in likelihood. Something is wrong");
      }
    }
    
    last_lb   = lb;
    last_like = like;
  }
  
  // Assertions
  if(arma::any(arma::vectorise(alpha_vb) < 0)) {stop("Negative values in alpha_vb");}
  if(arma::any(arma::vectorise(beta_vb) < 0)) {stop("Negative values in beta_vb");}
  
  // Compute and return the expected factors
  // E_W: Expectation of a Multinomial
  // E_H: Expectation of a Beta
  arma::mat E_W(F,Kmax);
  arma::mat E_H(N,Kmax);
  for(int f=0; f<F; f++){
    E_W.row(f) = gamma_vb.row(f) / arma::accu(gamma_vb.row(f)); 
  }
  for(int n=0; n<N; n++){
    E_H.row(n) = alpha_vb.row(n)  / (alpha_vb.row(n) + beta_vb.row(n));
  }
  
  // Assertions
  if(arma::any(arma::vectorise(E_W) < 0)) {stop("Negative values in E_W");}
  if(arma::any(arma::vectorise(E_H) < 0)) {stop("Negative values in E_H");}
  
  return Rcpp::List::create(Rcpp::Named("E_Z") = wrap(E_Z),
                            Rcpp::Named("E_W") = wrap(E_W),
                            Rcpp::Named("E_H") = wrap(E_H),
                            Rcpp::Named("lowerbounds") = wrap(lowerbounds),
                            Rcpp::Named("likelihoods") = wrap(likelihoods),
                            Rcpp::Named("kldivergences") = wrap(kldivergences));
}

