algo_VB_aspect <-function(V, Z, K,
                     alpha = 1, beta =1, gamma=1, 
                     iter=100){
  F <- dim(V)[1]
  N <- dim(V)[2]
  lowerbounds  <- rep(NA, iter)
  likelihoods <- rep(NA, iter)
  kldivergences <- rep(NA, iter)
  
  last_lowerbound <- -Inf
  
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
  cat("\nE_Z initialized")
  E_log_H     <- array(0, dim=c(N,K))
  E_log_1_H <- array(0, dim=c(N,K))
  E_log_W     <- array(0, dim=c(F,K))
  
  for(j in 1:iter){
    if(!j %% 20){cat('\n iteration:', j)}
    cat("\n*******************************")
    
    # Compute sums of Z marginals needed below
    E_Z_over_f <- t(colSums(E_Z, dims=1, na.rm = TRUE))
    E_Z_over_n <- rowSums(E_Z, dims=2, na.rm = TRUE) 
    
    masked_E_Z_plus  <- sweep(E_Z, c(1,3), V, '*', check.margin=FALSE)
    masked_E_Z_minus <- sweep(E_Z, c(1,3), 1-V, '*', check.margin=FALSE)
    E_Z_over_f_plus  <- t(colSums(masked_E_Z_plus,  dims=1, na.rm = TRUE))
    E_Z_over_f_minus <- t(colSums(masked_E_Z_minus, dims=1, na.rm = TRUE))
    
    
    gamma_vb <- gamma + E_Z_over_n       # F x K
    alpha_vb <- alpha + E_Z_over_f_plus  # N x K
    beta_vb  <- beta  + E_Z_over_f_minus # N x K
    
    # q_W ===================================================
    for(f in 1:F){
      E_log_W[f,] <- digamma(gamma_vb[f,]) - digamma(sum(gamma_vb[f,]))
    }
    if(j>1){
      lower <- lower_bound_aspect(E_Z, 
                                  E_log_W,
                                  E_log_H,
                                  E_log_1_H,
                                  gamma_vb, 
                                  alpha_vb, beta_vb, 
                                  alpha, beta, gamma, 
                                  V)
      cat(lower)
      if(lower < last_lowerbound) {
        stop("q(W): Decrease detected in lower bound. Something is wrong.")
      }
    }
    # q_H ===================================================
    for (n in 1:N){
      E_log_H[n,]   <- digamma(alpha_vb[n,]) - digamma(alpha_vb[n,] + beta_vb[n,])
      E_log_1_H[n,] <- digamma(beta_vb[n,]) -  digamma(alpha_vb[n,] + beta_vb[n,])
    }
    if(j>1){
      lower <- lower_bound_aspect(E_Z, 
                                  E_log_W,
                                  E_log_H,
                                  E_log_1_H,
                                  gamma_vb, 
                                  alpha_vb, beta_vb, 
                                  alpha, beta, gamma, 
                                  V)
      cat(lower)
      if(lower < last_lowerbound) {
        stop("q(H): Decrease detected in lower bound. Something is wrong.")
      }
    }
    
    # q_Z ===================================================
    for (n in 1:N){
      for(f in 1:F){
        if(is.na(V[f,n])){
          E_Z[f,,n] <- NA
        }
        
        E_Z[f,,n] <- exp(E_log_W[f,] + 
                           (V[f,n] * E_log_H[n,] + (1-V[f,n]) * E_log_1_H[n,]))
        E_Z[f,,n] <- E_Z[f,,n] / sum(E_Z[f,,n])
      }
    }
    
    if(j>1){
      lower <- lower_bound_aspect(E_Z, 
                                  E_log_W,
                                  E_log_H,
                                  E_log_1_H,
                                  gamma_vb, 
                                  alpha_vb, beta_vb, 
                                  alpha, beta, gamma, 
                                  V)
      cat(lower)
      if(lower < last_lowerbound) {
        stop("q(Z): Decrease detected in lower bound. Something is wrong.")
      }
    }
    
    # Compute lower bound ------------------------------------------------------
    if(TRUE){
      
      lowerR <- lower_bound_aspect_R(E_Z, 
                                    E_log_W,
                                    E_log_H,
                                    E_log_1_H,
                                    gamma_vb, 
                                    alpha_vb, beta_vb, 
                                    alpha, beta, gamma, 
                                    V)
      
      lower <- lower_bound_aspect(E_Z, 
                                   E_log_W,
                                   E_log_H,
                                   E_log_1_H,
                                   gamma_vb, 
                                   alpha_vb, beta_vb, 
                                   alpha, beta, gamma, 
                                   V)
      
      like <- likelihood_aspect(V,
                                  gamma_vb, 
                                  alpha_vb, beta_vb)
      
      cat("\n **Lower bound (R and Cpp):", lowerR, lower)
      cat("\n **Likelihood (Cpp):", like)
      cat("\n **KL (Cpp):", like - lower)
      lower <- lowerR
      if(lower < last_lowerbound) {
        stop("Decrease detected in lower bound. Somthing is wrong.")
      }
      lowerbounds[j] <- lower
      likelihoods[j] <- like
      kldivergences[j] <- like - lower
      last_lowerbound <- lower
    }
    
    # Store sample
    par(mfrow = c(2,1))
    plot(sort(apply(E_Z, 2, sum, na.rm=TRUE)), main=j, ylab = "Z column norms")
    plot(1:iter, lowerbounds, main=j, ylab = "lower bound")
    
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
              "E_H" = E_H,
              "lowerbounds" = lowerbounds,
              "likelihoods" = likelihoods,
              "kldivergences" = kldivergences)
  res
}

#################################################################
# Functions to monitor the the lower bound
#################################################################
E_log_p_W <- function(E_log_W, gamma){
  F <- nrow(E_log_W)
  K <- ncol(E_log_W)
  gamma <- rep(gamma, K)
  #x <- F*(lgamma(K*gamma) - K*lgamma(gamma))
  x <- 0
  for(f in 1:F){
    x <- x + lgamma(sum(gamma)) - sum(lgamma(gamma)) +
             sum((gamma-1) * E_log_W[f,])  
  }
  return(x)
}

E_log_q_W <- function(E_log_W, gamma_vb){
  F <- nrow(E_log_W)
  K <- ncol(E_log_W)
  x <- 0
  for(f in 1:F){
    x <- x + lgamma(sum(gamma_vb[f,])) - sum(lgamma(gamma_vb[f,])) +
             sum((gamma_vb[f,]-1)*E_log_W[f,])
  }
  return(x)
}

E_log_p_Z <- function(E_Z, E_log_W){
  sum(sweep(E_Z, c(1,2), E_log_W, '*', check.margin=FALSE), na.rm = TRUE)
}

E_log_q_Z <- function(E_Z){
  sum(E_Z * log(E_Z), na.rm = TRUE)
}

E_log_p_H <- function(E_log_H, E_log_1_H, alpha, beta){
  N <- nrow(E_log_H)
  K <- ncol(E_log_H)
  x <- N*K*(lgamma(alpha+beta) - lgamma(alpha) - lgamma(beta)) +
        (alpha-1)*sum(E_log_H) + (beta-1)*sum(E_log_1_H)
  return(x)
}

E_log_q_H <- function(E_log_H, E_log_1_H, alpha_vb, beta_vb){
  x <- sum(
    lgamma(alpha_vb + beta_vb) - lgamma(alpha_vb) - lgamma(beta_vb) +
    (alpha_vb-1) * E_log_H + (beta_vb-1) * E_log_1_H)
  return(x)
}

E_log_q_H_slow <- function(E_log_H, E_log_1_H, alpha_vb, beta_vb){
  N = nrow(E_log_H)
  K = ncol(E_log_H)
  x <- 0
  for(n in 1:N){
    for(k in 1:K){
      x <- x + 
        lgamma(alpha_vb[n,k] + beta_vb[n,k]) - 
          lgamma(alpha_vb[n,k]) - lgamma(beta_vb[n,k]) +
          (alpha_vb[n,k]-1) * E_log_H[n,k] + 
          (beta_vb[n,k]-1) * E_log_1_H[n,k]    
    }
  }
  return(x)
}


E_log_p_V <- function(V, E_Z, E_log_H, E_log_1_H){
  F <- nrow(V)
  N <- ncol(V)
  K <- ncol(E_log_H)
  V.probs <- matrix(NA, F,N)
  x <- 0
  for(n in 1:N){
    for(f in 1:F){
        V.probs[f,n] <- sum(E_Z[f,,n] * (V[f,n] * E_log_H[n,] +
                                     (1-V[f,n]) * E_log_1_H[n,]), 
                            na.rm = TRUE)
    }
  }
  x <- sum(V.probs, na.rm = TRUE)
  return(x)
}

lower_bound_aspect_R <- function(E_Z, 
                                 E_log_W,
                                 E_log_H,
                                 E_log_1_H,
                                 gamma_vb, 
                                 alpha_vb, beta_vb, 
                                 alpha, beta, gamma, 
                                 V){  
  
  F <- dim(V)[1]
  N <- dim(V)[2]
  K <- dim(E_Z)[2]
  lower <- 0
  
  # Use _rcpp version for debugging
  E_log_p_W <- E_log_p_W(E_log_W, gamma)  
  E_log_p_W_rcpp <- E_log_p_W_Rcpp(E_log_W, gamma)  
  
  E_log_q_W      <- E_log_q_W(E_log_W, gamma_vb)  
  E_log_q_W_rcpp <- E_log_q_W_Rcpp(E_log_W, gamma_vb)  
  
  E_log_p_Z <- E_log_p_Z(E_Z, E_log_W)
  E_log_p_Z_rcpp <- E_log_p_Z_Rcpp(E_Z, E_log_W, V)
  
  E_log_q_Z <- E_log_q_Z(E_Z)
  E_log_q_Z_rcpp <- E_log_q_Z_Rcpp(E_Z, V)
  
  E_log_p_H <- E_log_p_H(E_log_H, E_log_1_H, alpha, beta)
  E_log_p_H_rcpp <- E_log_p_H_Rcpp(E_log_H, E_log_1_H, alpha, beta)
  
  E_log_q_H <- E_log_q_H(E_log_H, E_log_1_H, alpha_vb, beta_vb)
  E_log_q_H_rcpp <- E_log_q_H_Rcpp(E_log_H, E_log_1_H, alpha_vb, beta_vb)
  
  E_log_p_V <- E_log_p_V(V, E_Z, E_log_H, E_log_1_H)
  E_log_p_V_rcpp <- E_log_p_V_Rcpp(V, E_Z, E_log_H, E_log_1_H)
  
  lower <- E_log_p_W + E_log_p_Z + E_log_p_H + E_log_p_V -
           E_log_q_W - E_log_q_Z - E_log_q_H
  cat("\n R lower bound", lower)
  # cat("\n R pW", E_log_p_W, E_log_p_W_rcpp)
  # cat("\n R qW", E_log_q_W, E_log_q_W_rcpp)
  # cat("\n R pH", E_log_p_H, E_log_p_H_rcpp)
  # cat("\n R qW", E_log_q_H, E_log_q_H_rcpp)
  # cat("\n R pZ", E_log_p_Z, E_log_p_Z_rcpp)
  # cat("\n R qZ", E_log_q_Z, E_log_q_Z_rcpp)
  # cat("\n R pV", E_log_p_V, E_log_p_V_rcpp)
  lower
}

