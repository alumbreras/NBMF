sample_gibbs_Z <- function(v_n, W, Z, alpha = 1, iter=100, burnin = 0.5){
  
  F <- dim(W)[1]
  K <- dim(W)[2]
  if(anyNA(W)) stop("W contains NAs")
  
  Z_samples <- array(NA, dim=c(nsamples, F, K))
  likelihood_matrix <- array(NA, dim=c(F, K))
  for(f in 1:F){
    likelihood_matrix[f,] <- W[f,]^(v_n[f]) * (1-W[f,])^(1-v_n[f])
  }

  for(j in 1:nsamples){
    for(f in 1:F){
      # Warning: we are also modeling the zeroes.
      # That's  why wz can't do thos
      # if(v_n[f]==0){
      #   Z[f,] <- 0
      # }
      Z[f,] <- 0 # remove current value
      L <- colSums(Z) # sum over f
      probs <- (alpha + L) * likelihood_matrix[f,] # compute probabilities
      probs <- probs/sum(probs)
      
      if(is.na(sum(probs))){
        cat('\n W_f', W[f,])
        cat('\n V_f', v_n[f])
        cat('\n likelihood_f', likelihood_matrix[f,])
      }
      k_chosen <- sample(K, size=1, prob=probs) # assign new component to f
      Z[f,k_chosen] <- 1 # re-assign value to chosen k
    }
    Z_samples[j,,] <- Z
  }
  #cat("\nZ_mean:",   colMeans(Z_samples, dims=1))
  colMeans(Z_samples, dims=1)
}


sample_gibbs_ZH <- function(v_n, W, Z, alpha = 1, iter=100, burnin = 0.5){
  
  F <- dim(W)[1]
  K <- dim(W)[2]
  if(anyNA(W)) stop("W contains NAs")
  
  Z_samples <- array(NA, dim=c(nsamples, F, K))
  likelihood_matrix <- array(NA, dim=c(F, K))
  for(f in 1:F){
    likelihood_matrix[f,] <- W[f,]^(v_n[f]) * (1-W[f,])^(1-v_n[f])
  }
  
  for(j in 1:nsamples){
    
    # Sample H
    alpha_post = rep(alpha, K) + colSums(Z)
    h = rdirichlet(1, alpha_post)
    
    # Sample Z
    for(f in 1:F){

      
      probs <- h * likelihood_matrix[f,] # computa probabilities
      probs <- probs/sum(probs)
      
      if(is.na(sum(probs))){
        cat('\n Warning!: SUM PROBS NULL!')
        cat('\n likelihood_f', likelihood_matrix[f,])
      }
      k_chosen <- sample(K, size=1, prob=probs) # assign new component to f
      Z[f,] <- 0 # remove current value
      Z[f,k_chosen] <- 1 # re-assign value to chosen k
    }
    Z_samples[j,,] <- Z
  }
  #cat("\nZ_mean:",   colMeans(Z_samples, dims=1))
  colMeans(Z_samples, dims=1)
}