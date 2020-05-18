Gibbs_DirBer <-function(V, Z, K=K, 
                                alpha=1, beta=1, gamma=1, 
                                iter=100, burnin = 1000){
  
  if(gamma < 1/K) {warning("Gamma very small")}
  
  # In Rcpp, the first label is 0
  Z <- Z - min(Z, na.rm=TRUE)
  
  res <- Gibbs_DirBer_finite_Rcpp(V, Z, K=K, 
                                    alpha=alpha, beta=beta, gamma=gamma, 
                                   iter=iter, burnin = burnin)
  E_Z <- res$E_Z
  cat("\nComputing E[W], E[H]")
  E_W <- compute_expectation_W_GibbsDirBer(res$Z_samples, gamma)
  E_H <- compute_expectation_H_GibbsDirBer(res$Z_samples, V, alpha, beta)
  Z_samples <- res$Z_samples
  
  res <- list(E_Z = E_Z,
              E_W = E_W,
              E_H = E_H,
              Z_samples = Z_samples,
              V = V,
              alpha=alpha,
              beta=beta,
              gamma=gamma)
  
  if(!all(res$Z_samples>=0, na.rm = TRUE)) {stop("E_Z contains negative values")}
  if(!all(res$E_Z>=0, na.rm = TRUE)) {stop("E_Z contains negative values")}
  if(!all(res$E_W>=0)) {stop("E_W contains negative values")}
  if(!all(res$E_H>=0)) {stop("E_H contains negative values")}
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
  V.train <- modelGibbsDirBer$V
  
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
    
    masked_Z_plus  <- sweep(Z, c(1,3), V.train, '*', check.margin=FALSE)
    masked_Z_minus <- sweep(Z, c(1,3), 1-V.train, '*', check.margin=FALSE)
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
  
  rownames(E_V) <- rownames(modelGibbsDirBer$V)
  colnames(E_V) <- colnames(modelGibbsDirBer$V)
  
  E_V
}


#' @title Draw a sample of W, H and V in Dirichlet Bernoulli model
#' @description Draw a sample of W, H and V=WH in Dirichlet Bernoulli model
#' given the samples of Z
#' @export
sample_VWH.GibbsDirBer <- function(modelGibbsDirBer){
  alpha <- modelGibbsDirBer$alpha
  beta  <- modelGibbsDirBer$beta
  gamma <- modelGibbsDirBer$gamma
  Z_samples <- modelGibbsDirBer$Z_samples
  
  # init from here to copy the row and col names
  V <- modelGibbsDirBer$V
  W <- modelGibbsDirBer$E_W
  H <- modelGibbsDirBer$E_H
  
  F <- dim(Z_samples)[1]
  N <- dim(Z_samples)[2]
  J <- dim(Z_samples)[3]
  K <- max(Z_samples, na.rm = TRUE)+1
  #W <- array(NA, dim=c(F,K))
  #H <- array(NA, dim=c(N,K))
  j <- sample((J-50):J,1) # use one of the last samples from Z
  
  cat("\n j:", j)
  Z <- Z_samples[,,j]
  Z <- matrix_to_tensor_R(Z, K)
  
  Z_over_f <- t(colSums(Z, dims=1, na.rm = TRUE)) # N x K
  Z_over_n <- rowSums(Z, dims=2, na.rm = TRUE)    # F x K
  
  masked_Z_plus  <- sweep(Z, c(1,3), V, '*', check.margin=FALSE)
  masked_Z_minus <- sweep(Z, c(1,3), 1-V, '*', check.margin=FALSE)
  Z_over_f_plus  <- t(colSums(masked_Z_plus,  dims=1, na.rm = TRUE))
  Z_over_f_minus <- t(colSums(masked_Z_minus, dims=1, na.rm = TRUE)) 
  
  # Sample W from Dirichlet
  for(f in 1:F){
    gamma_post <- gamma + Z_over_n[f,]
    gamma_post <- gamma_post/sum(gamma_post)
    W[f,] <- rdirichlet(1, gamma_post)
  }
  
  # Sample H from Beta
  for(n in 1:N){
    alpha_post <- alpha + Z_over_f_plus[n,]
    beta_post  <- beta + Z_over_f_minus[n,]
    H[n,] <- rbeta(K, alpha_post, beta_post)
  }
  
  # Sample V from Bernoulli(WH)
  probs <- W %*% t(H)
  for(f in 1:F){
    V[f,] <- rbinom(N, 1, prob = probs[f,])
  }
  
  list(W = W, H = H, V = V)
}

