#include <RcppArmadillo.h>
#include "commons.h"

using namespace Rcpp;
using namespace arma;


////////////////////////////////////////////////////////////////////////////////
// Counters for Z and C when there is one assignment per (f,n)
////////////////////////////////////////////////////////////////////////////////
struct counts {
  arma::imat f_ones; //declared as double to be ready for 'double' probabilities
  arma::imat f_zeros;
  arma::imat f_total;
  arma::imat n_total;
  arma::ivec fn_total;
  int Kmax;
  
  explicit counts(int f, int k, int n):  f_ones(n,k),
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
void compute_counts(const arma::imat& V, const arma::imat& Z, counts& cnt){
  
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
        int k = Z(f,n);
        cnt.f_ones(n,k)      += V(f,n);
        cnt.f_zeros(n,k)     += 1-V(f,n);
        cnt.f_total(n,k)     += 1;
        cnt.n_total(f,k)     += 1;
        cnt.fn_total(k)      += 1;
      }
    }
  }
}

// Assign k to a (f,n) pair and update counts
void set_assignment(int f, int n, int k, 
                    arma::imat& Z,
                    const arma::imat& V, counts& cnt){
  Z(f,n) = k;     
  cnt.f_ones(n,k)    += V(f,n);
  cnt.f_zeros(n,k)   += 1-V(f,n);
  cnt.f_total(n,k)   += 1;
  cnt.n_total(f,k)   += 1;
  cnt.fn_total(k)    += 1;
}

// Remove assigment from (f,n) pair and update counts
void remove_assignment(int f, int n, 
                       arma::imat& Z,
                       const arma::imat& V, counts& cnt){
  int k = Z(f,n);  // current assignment
  Z(f,n) = -1; // just as a safety measure
  cnt.f_ones(n,k)      -= V(f,n);
  cnt.f_zeros(n,k)     -= 1-V(f,n);
  cnt.f_total(n,k)     -= 1;
  cnt.n_total(f,k)     -= 1;
  cnt.fn_total(k)      -= 1;
}

////////////////////////////////////////////////////////////////////////////////
// Counters for Z and C when there are multiple assignments per (f,n)
////////////////////////////////////////////////////////////////////////////////
struct counts_multiple {
  arma::imat f_total;
  arma::imat n_total;
  arma::ivec fn_total;
  int Kmax;
  
  explicit counts_multiple(int f, int k, int n):  
    f_total(n,k),
    n_total(f,k),
    fn_total(k),
    Kmax(k) {
    f_total.zeros();
    n_total.zeros();
    fn_total.zeros();
  }
};

void compute_counts_multiple(const arma::imat& V, const arma::icube& C, 
                             counts_multiple& cnt){
  
  cnt.f_total.zeros();
  cnt.n_total.zeros();
  cnt.fn_total.zeros();
  
  int N = V.n_cols;
  int F = V.n_rows;
  int K = C.n_slices;
  arma::irowvec actives(K);
  
  for(int n=0; n<N; n++){
    for(int f=0; f<F; f++){
      if(V(f,n) != NA_INTEGER){
        actives = C.slice(n).row(f);
        cnt.f_total.row(n)     += actives;
        cnt.n_total.row(f)     += actives;
        cnt.fn_total           += actives;
      }
    }
  }
}

void set_assignment_multiple(int f, int n, arma::irowvec actives,
                             arma::icube& C,
                             const arma::imat& V, counts_multiple& cnt){
  C.slice(n).row(f) = actives;     
  cnt.f_total.row(n)   += actives;
  cnt.n_total.row(f)   += actives;
  cnt.fn_total         += actives;
}

void remove_assignment_multiple(int f, int n, 
                                arma::icube& C,
                                const arma::imat& V, counts_multiple& cnt){
  arma::irowvec actives = C.slice(n).row(f);  // current assignment
  C.slice(n).row(f).zeros();
  cnt.f_total.row(n)   -= actives;
  cnt.n_total.row(f)   -= actives;
  cnt.fn_total         -= actives;
}

////////////////////////////////////////////////////////////////////////////////
// Common expectation functions
////////////////////////////////////////////////////////////////////////////////

// Given a set of J matrix samples of Z (FxN)
// Return the expectation E[Z] where Z is a tensor representation
// with indicator vectors
arma::cube compute_expectation_Z(const arma::icube& Z_samples){
  int F = Z_samples.n_rows;
  int N = Z_samples.n_cols;
  int K = Z_samples.max();
  
  arma::icube Z_tensor(F,K+1,N);
  arma::cube E_Z(F,K+1, N);
  E_Z.zeros();
  
  int J = Z_samples.n_slices;
  for(int j=0; j<J; j++){
    matrix_to_tensor(Z_samples.slice(j), Z_tensor);
    E_Z += arma::conv_to<arma::cube>::from(Z_tensor);
  }
  E_Z = E_Z/J;
  return E_Z;
}

// Given a set of J matrix samples of Z (FxN)
// Return the expectation E[W | V]
// [[Rcpp::export]]
arma::mat compute_expectation_W_Rcpp(const arma::icube& Z_samples, double gamma){
  int F = Z_samples.n_rows;
  int K = Z_samples.max();
  arma::mat E_W(F,K+1);
  arma::cube E_Z = compute_expectation_Z(Z_samples);
  
  // Remove NA so that they are not counted in the sums
  E_Z.transform( [](double val) { return (std::isnan(val) ? double(0) : val);} );

  arma::mat E_Z_sufficient = sum(E_Z, 2); // sum over n
  
  for(int f=0; f<F; f++){
    E_W.row(f) = gamma + E_Z_sufficient.row(f);
  }
  for(int f=0; f<F; f++){
    E_W.row(f) /= arma::accu(E_W.row(f)); 
  }
  
  return E_W;
}
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////



// Given a set of J matrix samples of Z (FxN)
// Return the expectation E[W | V]
// [[Rcpp::export]]
arma::mat compute_expectation_H_dirdir_Rcpp(const arma::icube& C_samples, 
                                        double gamma){
  int N = C_samples.n_cols;
  int K = C_samples.max();
  arma::mat E_H(N,K+1);
  arma::cube E_C = compute_expectation_Z(C_samples);
  
  // Remove NA so that they are not counted in the sums
  E_C.transform( [](double val) { return (std::isnan(val) ? double(0) : val);} );
  arma::mat E_C_sufficient = sum(E_C, 0); // sum over f
  for(int n=0; n<N; n++){
    E_H.row(n) = gamma + E_C_sufficient.col(n).t();
  }
  for(int n=0; n<N; n++){
    E_H.row(n) /= arma::accu(E_H.row(n)); 
  }
  
  return E_H;
}


// Given a set of J matrix samples of Z (FxN) and V
// Return the expectation E[H | V]
// [[Rcpp::export]]
arma::mat compute_expectation_H_Rcpp(const arma::icube& Z_samples, 
                                const arma::imat& V,
                                double alpha, double beta){
  int J = Z_samples.n_slices;
  int F = Z_samples.n_rows;
  int N = Z_samples.n_cols;
  int K = Z_samples.max();
  arma::mat E_H(N,K+1);
  E_H.zeros();
  //arma::cube E_Z = compute_expectation_Z(Z_samples);
  //arma::mat E_Z_sufficient = sum(E_Z, 2); // sum over n

  for(int k=0; k<K; k++){
    for(int n=0; n<N; n++){
        E_H(n,k) = 0; 
      
        // Monte carlo approximation
        for(int j=0; j<J; j++){

          // Compute the counters over f
          int f_ones = 0;
          int f_total = 0;
          for(int f=0; f<F; f++){
            if(V(f,n) != NA_INTEGER){
              f_ones += (Z_samples.slice(j)(f,n)==k)*V(f,n);
              f_total += (Z_samples.slice(j)(f,n)==k);
              //Rcout << "f_ones: " <<  f_ones << std::endl;
              //Rcout << "f_total: " << f_total << std::endl;
            }
          }
          
          E_H(n,k) += (alpha + static_cast<double>(f_ones))/
                      (alpha + beta + static_cast<double>(f_total));
        }
        //Rcout << "E_H: " << n << ',' << k <<  " = "<< E_H(n,k) << std::endl;
        E_H(n,k) = E_H(n,k)/J; // mean
        //Rcout << "E_H/J: " << n << ',' << k <<  " = "<< E_H(n,k) << std::endl;
        
      }
  }
  // E_H = E_H/J;
  return E_H;
}

// [[Rcpp::export]]
SEXP Gibbs_DirBer_DP_Rcpp(const arma::imat& V, 
                       arma::imat Z,
                       double alpha = 1, double beta =1, double gamma=1, 
                       int iter=100,
                       double burnin = 0.5){

  double epsilon = 0.00000001;
  int Kmax = 200;
  int F = V.n_rows;
  int N = V.n_cols;
  double prob_nonempty;
  arma::rowvec	probs(Kmax);
  arma::rowvec	probs_BP(Kmax);
  arma::irowvec ans(Kmax);
  arma::rowvec	cumprobs(Kmax);
  arma::rowvec alpha_post;
  arma::rowvec beta_post;
  int nburnin = std::floor(iter*burnin);
  int nsamples = iter - std::floor(iter*burnin);
  arma::icube Z_samples(F,N,nsamples);
  Z_samples.zeros();
  
  int j = 0;

  // First counts
  // While sampling we only update individual elements
  counts cnt(F, Kmax, N);
  compute_counts(V,Z,cnt);

  // Resample the topic assigment k of each count
  for(int i=1; i<iter; i++){
    Rcout << "sample: " << i << '/' << iter << std::endl;
    
    for(int n=0; n<N; n++){
      for(int f=0; f<F; f++){

        // Remove current value
        remove_assignment(f, n, Z, V, cnt);
        
        // Compute probability of every component
        // the epsilons garantee that if counts f_ones f_zeros
        // go to zero, we don't get 0/0 but the correct value 0.5
        // (assuming alpha=beta)
        // if n_total.row(f) = 0, probability would be 0
        // but they will have a chance thanks to the CRP probability
        // which is proportional to gamma*0.5
        alpha_post = conv_to<rowvec>::from(cnt.f_ones.row(n)) + epsilon;
        beta_post  = conv_to<rowvec>::from(cnt.f_zeros.row(n)) + epsilon;
        probs_BP = (V(f,n) * alpha_post + (1- V(f,n)) * beta_post) / 
                   (alpha_post + beta_post);
        probs = conv_to<rowvec>::from(cnt.n_total.row(f)) % probs_BP;
        
        // Put probability gamma*0.5 in the first empty component
        // making it the representant of all empty components
        int k_free = std::distance(cnt.n_total.row(f).begin(), 
                                   std::find_if(cnt.n_total.row(f).begin(), 
                                                cnt.n_total.row(f).end(), 
                                                [](int x){ return x == 0;}));
        probs[k_free] = gamma*0.5;
        
        // Normalize and compute probabilities as usual
        probs = probs/accu(probs);
        
        // Draw assignment
        ans = ans.zeros();
        rmultinom(1, probs.begin(), Kmax, ans.begin());
        int k = arma::index_max(ans);
        
        set_assignment(f, n, k, Z, V, cnt);

      }
    }
    // After burning, store samples
    // (or sufficient statistics, if samples won't be needed)
    if(i >= nburnin){
      Z_samples.slice(j) = Z;
      j = j + 1;
    }
  }
  
  // Compute E[Z]
  arma::cube E_Z = compute_expectation_Z(Z_samples);
  
  return Rcpp::List::create(Rcpp::Named("Z_samples") = wrap(Z_samples),
                            Rcpp::Named("E_Z") = wrap(E_Z));
}


// [[Rcpp::export]]
SEXP Gibbs_DirBer_finite_Rcpp(const arma::imat& V, 
                              arma::imat Z, int K,
                              double alpha = 1, double beta =1, double gamma=1, 
                              int iter=100,
                              double burnin = 0.5){
  
  int Kmax = K;
  int F = V.n_rows;
  int N = V.n_cols;
  arma::vec	probs(Kmax);
  arma::vec	gamma_post(Kmax);
  arma::ivec ans(Kmax);
  arma::vec alpha_post(Kmax);
  arma::vec	beta_post(Kmax);
  int nburnin = std::floor(iter*burnin);
  int nsamples = iter - std::floor(iter*burnin);
  
  //Results storage
  arma::icube Z_samples(F,N,nsamples);
  Z_samples.zeros();

  int j = 0;
  
  // First counts
  // While sampling we only update individual elements
  counts cnt(F, Kmax, N);
  compute_counts(V,Z,cnt);
  
  // Resample the topic assigment k of each count
  for(int i=1; i<iter; i++){
    Rcout << "sample: " << i << '/' << iter << std::endl;
    
    for(int n=0; n<N; n++){
      for(int f=0; f<F; f++){
        if(V(f,n) != NA_INTEGER){
          
          remove_assignment(f, n, Z, V, cnt);
          
          // Compute probability of every component
          gamma_post =  gamma + conv_to<vec>::from(cnt.n_total.row(f));
          alpha_post = alpha + conv_to<vec>::from(cnt.f_ones.row(n));
          beta_post  = beta  + conv_to<vec>::from(cnt.f_zeros.row(n));
          probs = gamma_post % 
                      (V(f,n) * alpha_post + (1 - V(f,n)) * beta_post)/
                      (alpha_post + beta_post);
          probs = probs/accu(probs);
        
          // Draw assignment
          ans = ans.zeros();
          rmultinom(1, probs.begin(), Kmax, ans.begin());
          int k = arma::index_max(ans);
          
          set_assignment(f, n, k, Z, V, cnt);
        }
      }
    }
    
    // After burning, store samples
    // (or sufficient statistics, if samples won't be needed)
    //Rcout << "fn_total:" << cnt.fn_total.t() << std::endl;
    if(i >= nburnin){
      Z_samples.slice(j) = Z;
      j = j + 1;
    }
  }
  
  if(arma::any(arma::vectorise(alpha_post) < 0)) {stop("Negative values in alpha_post");}
  if(arma::any(arma::vectorise(beta_post) < 0)) {stop("Negative values in beta_post");}
  
  // Compute E[Z]
  arma::cube E_Z = compute_expectation_Z(Z_samples);
  
  return Rcpp::List::create(Rcpp::Named("Z_samples") = wrap(Z_samples),
                            Rcpp::Named("E_Z") = wrap(E_Z));
}


// [[Rcpp::export]]
SEXP Gibbs_DirDir_finite_Rcpp(const arma::imat& V, 
                                      arma::imat Z, int K,
                                      double alpha = 1, double gamma=1, 
                                      int iter=100,
                                      double burnin = 0.5){
   
   int Kmax = K;
   int F = V.n_rows;
   int N = V.n_cols;
   arma::vec	probs(Kmax);
   arma::ivec ans(Kmax);
   arma::vec gamma_post;
   arma::vec	alpha_post;
   int nburnin = std::floor(iter*burnin);
   int nsamples = iter - std::floor(iter*burnin);
   
   //Results storage
   arma::icube Z_samples(F,N,nsamples);
   arma::icube C_samples(F,N,nsamples);
   Z_samples.zeros();
   C_samples.zeros();
   
   arma::imat C = Z;
   
   int j = 0;
   
   // First counts
   // While sampling we only update individual elements
   counts cnt_z(F, Kmax, N);
   counts cnt_c(F, Kmax, N);
   compute_counts(V,Z,cnt_z);
   compute_counts(V,C,cnt_c);
   
   // Resample the topic assigment k of each count
   for(int i=1; i<iter; i++){
     Rcout << "sample: " << i << '/' << iter << std::endl;
     
     for(int n=0; n<N; n++){
       for(int f=0; f<F; f++){
         if(V(f,n) != NA_INTEGER){
           
           // If V_fn = 1, we sample Z and C together
           // Else, sample Z and then C
           if(V(f,n)){
             // Remove current assignment of Z(f,n), C(f,n)
             remove_assignment(f, n, Z, V, cnt_z);
             remove_assignment(f, n, C, V, cnt_c);
             gamma_post = gamma + conv_to<vec>::from(cnt_z.n_total.row(f));
             alpha_post = alpha + conv_to<vec>::from(cnt_c.f_total.row(n));
             
             // Compute probability for each component k
             probs = gamma_post % alpha_post;
             probs = probs/accu(probs);
  
             
             // Choose a component
             ans = ans.zeros();
             rmultinom(1, probs.begin(), Kmax, ans.begin());
             int k = arma::index_max(ans);
             set_assignment(f, n, k, Z, V, cnt_z);
             set_assignment(f, n, k, C, V, cnt_c);
             
           } else {
             // Sample Z(f,n) ------------------------------------------------
             // Remove current assignment of Z(f,n)
             remove_assignment(f, n, Z, V, cnt_z);
             gamma_post = gamma + conv_to<vec>::from(cnt_z.n_total.row(f));
             
             // Compute probability for each component k in Z(f,n)
             probs = gamma_post;
             probs(C(f,n))=0; //forbid the current component of C(f,n)
             probs = probs/accu(probs);
  
             // Choose a component for Z(f,n)
             ans = ans.zeros();
             rmultinom(1, probs.begin(), Kmax, ans.begin());
             int k_z = arma::index_max(ans);
   
             // Sample C(f,n) ------------------------------------------------
             // Remove current assignment of C(f,n)
             remove_assignment(f, n, C, V, cnt_c);
             alpha_post = alpha + conv_to<vec>::from(cnt_c.f_total.row(n));
             
             // Compute probability for each component k in C(f,n)
             probs = alpha_post;
             probs(k_z)=0; //forbid the current component of Z(f,n)
             probs = probs/accu(probs);
  
             // Choose a component for C(f,n)
             ans = ans.zeros();
             rmultinom(1, probs.begin(), Kmax, ans.begin());
             int k_c = arma::index_max(ans);
  
             set_assignment(f, n, k_z, Z, V, cnt_z);
             set_assignment(f, n, k_c, C, V, cnt_c);
  
           }
         }
       }
     }
     // After burning, store samples
     // (or sufficient statistics, if samples won't be needed)
     //Rcout << "fn_total:" << cnt_z.fn_total.t() << std::endl;
     if(i >= nburnin){
       Z_samples.slice(j) = Z;
       C_samples.slice(j) = C;
       j = j + 1;
     }
   }
   
   // Compute E[Z]
   arma::cube E_Z = compute_expectation_Z(Z_samples);
   arma::cube E_C = compute_expectation_Z(C_samples);
   
   return Rcpp::List::create(Rcpp::Named("Z_samples") = wrap(Z_samples),
                             Rcpp::Named("C_samples") = wrap(C_samples),
                             Rcpp::Named("E_Z") = wrap(E_Z),
                             Rcpp::Named("E_C") = wrap(E_C)
                             );
}

// Take J samples from W, H and then V
// using the samples of Z
// [[Rcpp::export]]
SEXP samples_VWH_DirBer(const arma::imat& V, const arma::icube& Z_samples, 
                        double gamma, double alpha, double beta){

  int J = Z_samples.n_slices;
  int F = Z_samples.n_rows;
  int N = Z_samples.n_cols;
  int K = Z_samples.max()+1;
  
  arma::vec	gamma_post(K);
  arma::vec alpha_post(K);
  arma::vec	beta_post(K);
  arma::mat probs(F,N);
  
  arma::cube W_samples(F,K,J);
  arma::cube H_samples(N,K,J);
  arma::icube V_samples(F,N,J);
  
  counts cnt(F, K, N);
  
  for(int j=0; j<J; j++){
    Rcout << "sample: " << j << '/' << J << std::endl;
    
    compute_counts(V, Z_samples.slice(j), cnt);

    // Sample W from Dirichlet (using normalized Gammas)
    for(int f=0; f<F; f++){
      gamma_post = gamma + conv_to<vec>::from(cnt.n_total.row(f));
      for(int k=0; k<K; k++){
        W_samples.slice(j)(f,k) = arma::randg(1, arma::distr_param(gamma_post[k],1.0))[0];
      }
      W_samples.slice(j).row(f) /= arma::accu(W_samples.slice(j).row(f));
    }
    
    // Sample H from Beta
    for(int n=0; n<N; n++){
      alpha_post = alpha + conv_to<vec>::from(cnt.f_ones.row(n));
      beta_post  = beta  + conv_to<vec>::from(cnt.f_zeros.row(n));
      for(int k=0; k<K; k++){
        H_samples.slice(j)(n,k) = Rcpp::as<double>(Rcpp::rbeta(1, alpha_post[k], beta_post[k]));
      }
    }
    
    // Sample V from Bernoulli(WH)
    // TODO: gives wrong results
    probs = W_samples.slice(j) * H_samples.slice(j).t();
    for(int f=0; f<F; f++){
      for(int n=0; n<N; n++){
        if(Rcpp::as<double>(Rcpp::runif(1)) < probs(f,n)){
          V_samples.slice(j)(f,n) = 1;
        } else{
          V_samples.slice(j)(f,n) = 0;
        }
      }
    }
      
  }
  
  return Rcpp::List::create(Rcpp::Named("W_samples") = wrap(W_samples),
                            Rcpp::Named("H_samples") = wrap(H_samples),
                            Rcpp::Named("V_samples") = wrap(V_samples));
}