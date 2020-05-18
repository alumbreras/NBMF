# This implements MCEM inference for a Dirichlet Bernouilli NMF

#' @title NMF with MMLE
#' @param V matrix to be factorized
#' @param K number of latent factors (or dimensions)
#' @param W initial matrix W
#' @param maxiters maximum number of iterations
#' @param nsamples number of samples per iteration in the E-step
#' @param mc mcem or saem.
#' @details The number of iterations is set to maxiters.
#' In the future, there should be a convergence check in case
#' the KL converges before maxiters.
#' @return W last estimator of W
#'         W_traces a K x samples - matrix vector with the evolution of row norms 
mmle_dirber_mcem <- function(V, initW, alpha=1, 
                          maxiters = 100, nsamples = 50, burnin=0.5){
  
  F <- dim(V)[1]
  K <- dim(initW)[2]
  N <- dim(V)[2]
  npostburnin <- floor(nsamples*burnin) # burnin samples
  
  if(dim(V)[1] != dim(initW)[1]) {
    stop("Dimensions between V and initW do not match")
  }
  
  if(nsamples<2) stop('Too few samples by iteration')
  
  C_samples_mean <- array(0, dim=c(F, K, N))
  W_traces <- array(NA, dim=c(maxiters, F, K))
  times <- rep(NA, maxiters)
  
  W <- initW
  
  # Initialize the first latent matrix with everything in the first dimension
  C_samples_mean[,1,] = V
  #C_samples_mean[,1,] = 1
  
  # Initialize the first matrix using row-wise multinomials
  # so that we have a valid matrix where sum_k C = V
  # TODO: slow in big data. 
  # for(f in 1:F){
  #   cat("\n f:", f)
  #   for(n in 1:N){
  #     #C_samples_mean[f,,n] <- t(rmultinom(1, size=V[f,n], prob=rep(1/K, K)))[1,]
  #     C_samples_mean[f,,n] <- V[f,n]/K
  #   }
  # }
  
  niter <- 1
  idx_trace <- 1
  for (niter in 1:maxiters){
    
    cat("\n ***************")
    cat("\n * niter",niter, "/", maxiters)
    cat("\n ***************")
    
    # Gibbs --------------------------------------------------------------------
    
    cat("\n E-step")
    # E-step
    for(n in 1:N){
      if(!n %% 20){cat('\n iter:', niter, '. Gibbs in column ', n, "/", N)}
      samples <- sample_gibbs_Z(V[,n], W, C_samples_mean[,,n],
                                  alpha=alpha,
                                  iter=nsamples, burnin=0)
      C_samples_mean[,,n] <- samples
    }
    
    
    # M-step
    #TODO check that components with 0 assignments obtain \sum W_f=0
    E_Z <- C_samples_mean
    for(f in 1:F){
      for(k in 1:K){
        #cat('\n', f,k, sum(E_Z[f,k,]*V[f,]))
        #cat('\n',sum(E_Z[f,k,]))
        W[f,k] <- sum(E_Z[f,k,]*V[f,])/sum(E_Z[f,k,])
        #cat('\n', f,k,'.........', W[f,k])
      }
    }
    
    # If some feature-topic is empty for all n,
    # W[f,k] is 0.
    # Since sum(E_Z[f,k,])=0, the above division gives 0/0 = NA
    W[is.na(W)] <- 0
    
    if(anyNA(W)){
      stop("W contains NAs")
    }
    
    cat("\ncolSums W\n:", colSums(W))
    cat("\ncolSumns Z\n:", apply(C_samples_mean, 2, function(x) sum(x)))
    
    W_traces[niter,,] <- W
    times[niter] <- Sys.time()
    
    
  } # end iters
  
  times <- times - min(times)
  list(V=V, W=W, W_traces = W_traces, times = times)
}