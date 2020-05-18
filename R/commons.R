#' @export
loglikelihood <- function(x, ...){ 
  UseMethod("loglikelihood", x)
}

#' @export
expectation <- function(x){ 
  UseMethod("expectation", x)
}

#' @title Dirichlet-Bernoulli loglikelihood
#' @description Compute likelihood of some matrix V  p(V* | V)
#' @param V A test matrix V. Same as the original V but with NA
#' in the training values
#' @details It uses the samples from WH and p(V* | V) = E[WH]
#' @export
loglikelihood.Bernoulli <- function(modelBernoulli, V_){
  E_V = expectation(modelBernoulli)
  V.probs <- V_*E_V + (1-V_)*(1-E_V)
  like <- sum(log(V.probs), na.rm = TRUE)
  res <- list(loglikelihood = like,
              V.probs = V.probs)
  return(res)
}


#' @title Expectation of Variational Inference model
#' @param E_W F x K matrix with expected dictionary W
#' @param E_H N x K matrix with expected activations H
#' @export
expectation.VBBernoulli <- function(modelVB){
  modelVB$E_W %*% t(modelVB$E_H) 
}
