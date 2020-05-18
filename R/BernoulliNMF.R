#' @title Collapsed Gibbs Sampler for Bernoulli models
#' @param V Observed matrix with FxN dimensions
#' @param K Number of latent dimensions
#' @param alpha parameter for the Beta distribution
#' @param beta parameter for the Beta distribution
#' @param gamma parameter for the Dirichlet distribution
#' @param iter number of Gibbs samples
#' @param burnin proportions of initial samples that are ignored for burnin
#' @description Samples the auxiliary variables (Z or Z,C) using a collapsed
#' Gibbs sampler.
#' @return FxNxJ arrays with samples from the variables
BernoulliNMF <- function(V, K=K, K_init = K,
                              model="DirBer", method="Gibbs",
                              alpha=1, beta=1, gamma=1, 
                              iter=1000, burnin = 0.5,
                              Z_init = NULL){
  
  if(method == "Gibbs"){
    K_init <- 10
  }
  
  # Random initialization of Z matrix
  if(is.null(Z_init)){
    F <- dim(V)[1]
    N <- dim(V)[2]
    Z <- array(NA, dim=c(F, N))
    for(n in 1:N){
      for(f in 1:F){
        if(!is.na(V[f,n])){
          Z[f,n] <- sample(K_init,1)
        }
      }
    }
  } else{
    stopifnot(dim(Z)[1] ==  dim(V)[1])
    stopifnot(dim(Z)[2] ==  dim(V)[2])
  }


  # Gibbs samplers -------------------------------------------------------------
  if(method == "Gibbs"){
    
    # Compute some iterations of VB and use E_Z as the final initialization for Z
    res <- VB_aspectRcpp(V, Z, K=K_init, alpha = alpha, beta=beta, gamma=gamma, iter=100)
    E_Z <- res$E_Z
    for(n in 1:N){
      for(f in 1:F){
        if(!is.na(V[f,n])){
          Z[f,n] <- which.max(E_Z[f,,n])
        }
      }
    }
    
    if(model == "DirBer"){
      res <- Gibbs_DirBer(V, Z, K=K, 
                         alpha=alpha, beta=beta, gamma=gamma, 
                         iter=iter, burnin = burnin)
    }
    
    if(model == "DirBerDP"){
      res <- Gibbs_DirBer_DP(V, Z, 
                          alpha=alpha, beta=beta, gamma=gamma, 
                          iter=iter, burnin = burnin)
    }
    
    if(model == "DirDir"){
      res <- Gibbs_DirDir(V, Z, K=K, 
                          alpha=alpha, gamma=gamma, 
                          iter=iter, burnin = burnin)
    }
  }
  
  # Variational Inference ------------------------------------------------------
  if(method == "VB"){
    
    if(model == "aspect"){
      # VB as in the Aspect Model paper (i.e. uncollapsed)
      cat("\nRunning VB_aspect")
      res <- VB_aspect(V, Z, K, alpha = alpha, beta=beta, gamma=gamma, iter=iter)
    }
    
    if(model == "aspectRcpp"){
      # VB as in the Aspect Model paper (i.e. uncollapsed)
      cat("\nRunning VB_aspect")
      res <- VB_aspectRcpp(V, Z, K, alpha = alpha, beta=beta, gamma=gamma, iter=iter)
    }
    
    if(model == "DirBer"){
      cat("\nRunning VB_DirBer")
      #stop('Collapsed version. Not implemented (buggy, indeed)')
      # Collapsed VB with truncated Dirichlet prior
      res <- VB_DirBer(V, Z, K, alpha=alpha, beta=beta, gamma=gamma, iter=iter)
    }
    
    if(model == "DirDir"){
      # Collapsed VB with truncated Dirichlet prior
      stop('Collapsed version. Not implemented')
    }
  }
  
  if(dim(res$E_W)[1] != dim(V)[1]) { stop("Dimensions E_W and V do not match") }
  if(dim(res$E_H)[1] != dim(V)[2]) { stop("Dimensions E_H and V do not match") }
  rownames(res$E_W) <- rownames(V)
  rownames(res$E_H) <- colnames(V)
  res
}
