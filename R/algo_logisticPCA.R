# WARNING: competitive algos to optimize the logisticPCA are not that simple
# I would better use libraries and forget about that.

#' @title Logistic PCA
#' @description Collapsed Variational Inference under the Truncated Symmetric Dirichlet Prior
#' Instead of computing the real limit K -> inf, we keep the origial Dirichlet-Multinomial
#' with a large K and Dirichlet parameters alpha/K
algo_logisticPCA <-function(V, K, iter=100){
  
  
  F <- dim(V)[1]
  N <- dim(V)[2]
  
  if(K < max(Z, na.rm = TRUE)){
    stop("K smaller than the initial number of components in Z")
  }
  
  # Z tensor represents E[Z], keeping the last E[Z[f,k,n]] of every element
  E_Z <- array(0, dim=c(F,K,N))
  for(f in 1:F){
    for(n in 1:N){
      k <- Z[f,n]
      if(!is.na(V[f,n])){
        E_Z[f,k,n] <- 1
      }
    }
  }
  
  
  # Matrices for bookkeeping of counts
  #L_plus    <- array(NA, dim=c(N,K))
  #L_minus   <- array(NA, dim=c(N,K))
  #Z_minus_n <- array(NA, dim=c(F,K))
  
  # Initial counts
  E_Z_over_f <- t(colSums(E_Z, dims=1, na.rm = TRUE)) # N x K
  E_Z_over_n <- rowSums(E_Z, dims=2, na.rm = TRUE)    # F x K
  masked_E_Z_plus  <- sweep(E_Z, c(1,3), V, '*', check.margin=FALSE)
  masked_E_Z_minus <- sweep(E_Z, c(1,3), 1-V, '*', check.margin=FALSE)
  E_Z_over_f_plus  <- t(colSums(masked_E_Z_plus,  dims=1, na.rm = TRUE))
  E_Z_over_f_minus <- t(colSums(masked_E_Z_minus, dims=1, na.rm = TRUE)) 
  
  j <- 1
  n <- 1
  f <- 1
  for(j in 1:iter){
    cat("\n j:", j)
    for (n in 1:N){
      for(f in 1:F){
        
        if(is.na(V[f,n])){
          next
        }
        
        # remove current expectation values from counters
        current_val    <- E_Z[f,,n]
        E_Z_over_f[n,]       <- E_Z_over_f[n,] - current_val
        E_Z_over_n[f,]       <- E_Z_over_n[f,] - current_val
        E_Z_over_f_plus[n,]  <- E_Z_over_f_plus[n,]  - current_val * V[f,n]
        E_Z_over_f_minus[n,] <- E_Z_over_f_minus[n,] - current_val * (1-V[f,n])
        E_Z[f,,n] <- 0
        
        gamma_vb <- gamma + E_Z_over_n       # F x K
        alpha_vb <- alpha + E_Z_over_f_plus  # N x K
        beta_vb  <- beta  + E_Z_over_f_minus # N x K
        
        probs <- gamma_vb[f,] * 
          (V[f,n] * alpha_vb[n,] + (1-V[f,n]) * beta_vb[n,]) /
          (alpha_vb[n,] + beta_vb[n,])
        
        
        probs <- probs / sum(probs)
        
        # Update expectation of each k in Z[f,,n], and counters
        E_Z[f,,n] <- probs 
        E_Z_over_f[n,]       <- E_Z_over_f[n,] + probs
        E_Z_over_n[f,]       <- E_Z_over_n[f,] + probs
        E_Z_over_f_plus[n,]  <- E_Z_over_f_plus[n,]  + probs * V[f,n]
        E_Z_over_f_minus[n,] <- E_Z_over_f_minus[n,] + probs * (1-V[f,n])
      }
    }
    plot(sort(colSums(rowSums(E_Z, dims = 2, na.rm = TRUE))))
  }
  
  cat("\n Computing expectations:")
  # Compute and return the expected factors
  # E_W: Expectation of a Multinomial
  # E_H: Expectation of a Beta
  E_H     <- array(0, dim=c(N,K))
  E_W     <- array(0, dim=c(F,K))
  for(f in 1:F){
    cat("\n Computing expectations f:", f)
    E_W[f,] <- gamma_vb[f,]/sum(gamma_vb[f,]) 
  }
  for(n in 1:N){
    cat("\n Computing expectations n:", n)
    E_H[n,] <- alpha_vb[n,]  / (alpha_vb[n,] + beta_vb[n,])
  }
  
  res <- list("E_Z" = E_Z,
              "E_W" = E_W,
              "E_H" = E_H)
  res
}
