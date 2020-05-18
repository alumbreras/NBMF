#' @title Collapsed VB for the Dirichlet-Bernoulli model
#' @description Collapsed Variational Inference under the Truncated Symmetric Dirichlet Prior
#' Instead of computing the real limit K -> inf, we keep the origial Dirichlet-Multinomial
#' with a large K and Dirichlet parameters alpha/K
VB_DirBer <-function(V, Z, K, alpha = 1, beta =1, gamma=1, iter=100){
  
  # In Rcpp, the first label is 0
  Z <- Z - min(Z, na.rm=TRUE)
  
  res <- VB_DirBer_Rcpp(V, Z, K, alpha=alpha, beta=beta, gamma=gamma, iter=iter)
  
  if(!all(res$E_W>=0)) {stop("E_W contains negative values")}
  if(!all(res$E_H>=0)) {stop("E_W contains negative values")}
  class(res) <- c("VBDirBer", "VBBernoulli", "Bernoulli")

  res
}
