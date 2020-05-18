Gibbs_DirBer_DP <-function(V, Z, 
                        alpha=1, beta=1, gamma=1, 
                        iter=100, burnin = 1000){
  
  # In Rcpp, the first label is 0
  Z <- Z - min(Z, na.rm=TRUE)
  
  res <- Gibbs_DirBer_DP_Rcpp(V, Z, 
                              alpha=alpha, beta=beta, gamma=gamma, 
                              iter=iter, burnin = burnin)
  
  E_Z <- res$E_Z
  E_W <- compute_expectation_W_GibbsDirBer(res$Z_samples, gamma)
  E_H <- compute_expectation_H_GibbsDirBer(res$Z_samples, V, alpha, beta)
  Z_samples <- res$Z_samples
  
  res <- list(E_Z = E_Z,
              E_W = E_W,
              E_H = E_H,
              Z_samples = Z_samples,
              alpha=alpha,
              beta=beta,
              gamma=gamma)
  if(!all(E_W>=0)) {stop("E_W contains negative values")}
  if(!all(E_H>=0)) {stop("E_W contains negative values")}
  class(res) <- c("GibbsDirBer", "Bernoulli")
  res
}

compute_expectation_W_GibbsDirBer <- function(Z_samples, gamma){
  E_W <- compute_expectation_W_Rcpp(Z_samples, gamma)
  E_W
}

compute_expectation_H_GibbsDirBer <- function(Z_samples, V, alpha, beta){
  E_H <- compute_expectation_H_Rcpp(Z_samples, V, alpha, beta)
  E_H
}

#' @title Expectation in Dirichlet Bernoulli model
#' @description Compute expected matrix V
#' @details It gets V=E[WH] given each sample from Z
#' @export
expectation.GibbsDirBer <- function(modelGibbsDirBer){
  alpha <- modelGibbsDirBer$alpha
  beta  <- modelGibbsDirBer$beta
  gamma <- modelGibbsDirBer$gamma
  Z_samples <- modelGibbsDirBer$Z_samples
  F <- dim(Z_samples)[1]
  N <- dim(Z_samples)[2]
  J <- dim(Z_samples)[3]
  K <- max(Z_samples, na.rm = TRUE)+1
  E_W <- array(NA, dim=c(F,K))
  E_H <- array(NA, dim=c(N,K))
  E_V <- array(0, dim=c(F,N))
  for(j in 1:J){
    cat("\n j:", j)
    Z <- Z_samples[,,j]
    Z <- matrix_to_tensor_R(Z, K)
    
    Z_over_f <- t(colSums(Z, dims=1, na.rm = TRUE)) # N x K
    Z_over_n <- rowSums(Z, dims=2, na.rm = TRUE)    # F x K
    
    masked_Z_plus  <- sweep(Z, c(1,3), V, '*', check.margin=FALSE)
    masked_Z_minus <- sweep(Z, c(1,3), 1-V, '*', check.margin=FALSE)
    Z_over_f_plus  <- t(colSums(masked_Z_plus,  dims=1, na.rm = TRUE))
    Z_over_f_minus <- t(colSums(masked_Z_minus, dims=1, na.rm = TRUE)) 
    
    # Expected W
    for(f in 1:F){
      E_W[f,] <- gamma + Z_over_n[f,]
      E_W[f,] <- E_W[f,]/sum(E_W[f,])
    }
    
    # Expected H
    for(n in 1:N){
      E_H[n,] <- (alpha + Z_over_f_plus[n,]) / (alpha + beta + Z_over_f[n,])
    }
    
    # Expected V
    E_V <- E_V + E_W %*% t(E_H)
  }
  E_V <- E_V/J
  E_V
}


