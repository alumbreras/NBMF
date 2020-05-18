Gibbs_DirDir <-function(V, Z, K=K, 
                                           alpha=1, gamma=1, 
                                           iter=100, burnin = 1000){
  # In Rcpp, the first label is 0
  Z <- Z - min(Z, na.rm=TRUE)
  
  res <- Gibbs_DirDir_finite_Rcpp(V, Z, K=K, 
                                         alpha=alpha, gamma=gamma, 
                                         iter=iter, burnin = burnin)
  E_Z <- res$E_Z
  E_C <- res$E_C
  cat("\nComputing E[W], E[H]")
  E_W <- compute_expectation_W_GibbsDirDir(res$Z_samples, gamma)
  E_H <- compute_expectation_H_GibbsDirDir(res$C_samples, alpha)
  Z_samples <- res$Z_samples
  C_samples <- res$C_samples

  res <- list(E_Z = E_Z,
              E_C = E_C,
              E_W = E_W,
              E_H = E_H,
              V = V,
              Z_samples = Z_samples,
              C_samples = C_samples,
              alpha=alpha,
              gamma=gamma)

  if(!all(E_W>=0)) {stop("E_W contains negative values")}
  if(!all(E_H>=0)) {stop("E_W contains negative values")}
  if(dim(E_W)[1] != dim(V)[1]) { stop("Dimensions E_H and V do not match") }
  if(dim(E_H)[1] != dim(V)[2]) { stop("Dimensions E_H and V do not match") }
  class(res) <- c("GibbsDirDir", "Bernoulli")
  res
}

compute_expectation_W_GibbsDirDir <- function(Z_samples, gamma){
  E_W <- compute_expectation_W_Rcpp(Z_samples, gamma)
  E_W
}

compute_expectation_H_GibbsDirDir <- function(C_samples, alpha){
  E_H <- compute_expectation_H_dirdir_Rcpp(C_samples, alpha)
  E_H
}


#' @export
expectation.GibbsDirDir <- function(modelGibbsDirDir){
  alpha <- modelGibbsDirDir$alpha
  gamma <- modelGibbsDirDir$gamma
  Z_samples <- modelGibbsDirDir$Z_samples
  C_samples <- modelGibbsDirDir$C_samples
  #V <- modelGibbsDirBer$V
  
  F <- dim(Z_samples)[1]
  N <- dim(Z_samples)[2]
  J <- dim(Z_samples)[3]
  
  # E_W, E_H set to the same max K. If a dimension is not used by either W or H
  # then its final contribution to V will be zero.
  Kz <- max(Z_samples, na.rm = TRUE)+1
  Kc <- max(C_samples, na.rm = TRUE)+1
  K <- max(Kz, Kc) 
  E_W <- array(NA, dim=c(F,K))
  E_H <- array(NA, dim=c(N,K))
  E_V <- array(0, dim=c(F,N))
  for(j in 1:J){
    cat("\n j:", j)
    Z <- Z_samples[,,j]
    C <- C_samples[,,j]
    Z <- matrix_to_tensor_R(Z, K)
    C <- matrix_to_tensor_R(C, K)
    
    Z_over_f <- t(colSums(Z, dims=1, na.rm = TRUE)) # N x K
    Z_over_n <- rowSums(Z, dims=2, na.rm = TRUE)    # F x K
    C_over_f <- t(colSums(C, dims=1, na.rm = TRUE)) # N x K
    C_over_n <- rowSums(C, dims=2, na.rm = TRUE)    # F x K
    
    # Expected W
    for(f in 1:F){
      E_W[f,] <- gamma + Z_over_n[f,]
      E_W[f,] <- E_W[f,]/sum(E_W[f,])
    }
    
    # Expected H
    for(n in 1:N){
      E_H[n,] <- alpha + C_over_f[n,]
      E_H[n,] <- E_H[n,]/sum(E_H[n,])
    }
    
    # Expected V
    E_V <- E_V + E_W %*% t(E_H)
  }
  E_V <- E_V/J
  rownames(E_V) <- rownames(modelGibbsDirDir$V)
  colnames(E_V) <- colnames(modelGibbsDirDir$V)
  
  E_V
}


#' @title Draw a sample of W, H and V in Dirichlet Dirichlet model
#' @description Draw a sample of W, H and V=WH in Dirichlet Bernoulli model
#' given the samples of Z
#' @export
sample_VWH.GibbsDirDir <- function(modelGibbsDirDir){
  alpha <- modelGibbsDirDir$alpha
  gamma <- modelGibbsDirDir$gamma
  Z_samples <- modelGibbsDirDir$Z_samples
  C_samples <- modelGibbsDirDir$C_samples

  F <- dim(Z_samples)[1]
  N <- dim(Z_samples)[2]
  J <- dim(Z_samples)[3]
  Kz <- max(Z_samples, na.rm = TRUE)+1
  Kc <- max(C_samples, na.rm = TRUE)+1
  K <- max(Kz, Kc) 
  
  # init from here to copy the row and col names
  V <- modelGibbsDirDir$V
  W <- modelGibbsDirDir$E_W
  H <- modelGibbsDirDir$E_H

  j <- sample((J-50):J,1) # use one of the last samples from Z
  
  cat("\n j:", j)
  Z <- Z_samples[,,j]
  C <- C_samples[,,j]
  Z <- matrix_to_tensor_R(Z, K)
  C <- matrix_to_tensor_R(C, K)
  
  Z_over_f <- t(colSums(Z, dims=1, na.rm = TRUE)) # N x K
  Z_over_n <- rowSums(Z, dims=2, na.rm = TRUE)    # F x K
  C_over_f <- t(colSums(C, dims=1, na.rm = TRUE)) # N x K
  C_over_n <- rowSums(C, dims=2, na.rm = TRUE)    # F x K

  # Sample W from Dirichlet
  for(f in 1:F){
    gamma_post <- gamma + Z_over_n[f,]
    gamma_post <- gamma_post/sum(gamma_post)
    W[f,] <- rdirichlet(1, gamma_post)
  }
  
  # Sample H from Beta
  for(n in 1:N){
    alpha_post <- alpha + C_over_f[n,]
    alpha_post <- alpha_post/sum(alpha_post)
    H[n,] <- rdirichlet(1, alpha_post)
  }
  
  # Sample V from Bernoulli(WH)
  probs <- W %*% t(H)
  for(f in 1:F){
    V[f,] <- rbinom(N, 1, prob = probs[f,])
  }
  
  list(W = W, H = H, V = V)
}


