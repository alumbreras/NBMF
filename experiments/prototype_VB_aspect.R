# Easy, synthetic example with a block diagonal matrix to test the 
# self-regularization property of the model

devtools::load_all()
library(dplyr)
library(tidyr)
library(ggplot2)
library(Matrix)
library(gtools)
library(data.table)

text_components <-  7
text_features <- 7

load("../data/adjacency_parlamentcat_sorted.RData")

V <- V[,-1]

K <- 50
F <- dim(V)[1]
N <- dim(V)[2]
Z <- array(0, dim=c(F, N))
for(n in 1:N){
  for(f in 1:F){
    Z[f,n] <- sample(K,1)
  }
}
Z_matrix <- Z


K = 50
alpha = 1
beta = 1
gamma = 1#/K

F = dim(Z_matrix)[1]
N = dim(Z_matrix)[2]
if(K < max(Z_matrix)){
  stop("K smaller than the initial number of components in Z")
}



# Z tensor represents E[Z], keeping the last E[Z[f,k,n]] of every element
E_Z <- array(0, dim=c(F,K,N))
for(f in 1:F){
  for(n in 1:N){
    k <- Z_matrix[f,n]
    E_Z[f,k,n] <- 1
  }
}

E_log_h     <- array(0, dim=c(K,N))
E_log_1_h_n <- array(0, dim=c(K,N))
E_log_w     <- array(0, dim=c(F,K))


E_Z <- VB_aspect(V, Z_matrix, K, alpha = 1, beta =1, gamma=1, iter=100)

j <- 1
n <- 1
f <- 1
iter <- 1000
for(j in 1:iter){
  if(!j %% 20){cat('\n iteration:', j)}
  cat("\n*******************************")
  
  # Compute sums of Z marginals needed below
  E_Z_over_f <- t(colSums(E_Z, dims=1))
  E_Z_over_n <- rowSums(E_Z, dims=2) 
  
  masked_E_Z_plus  <- sweep(E_Z, c(1,3), V, '*', check.margin=FALSE)
  masked_E_Z_minus <- sweep(E_Z, c(1,3), 1-V, '*', check.margin=FALSE)
  E_Z_over_f_plus  <- t(colSums(masked_E_Z_plus,  dims=1))
  E_Z_over_f_minus <- t(colSums(masked_E_Z_minus, dims=1))
  
  
  gamma_vb <- gamma + E_Z_over_n       # F x K
  alpha_vb <- alpha + E_Z_over_f_plus  # N x K
  beta_vb  <- beta  + E_Z_over_f_minus # N x K
  
  # q_W ===================================================
  for(f in 1:F){
    E_log_w[f,] <- digamma(gamma_vb[f,]) -
                   digamma(sum(gamma_vb[f,]))
  }
  
  # q_H ===================================================
  for (n in 1:N){
    E_log_h[,n]     <- digamma(alpha_vb[n,]) - digamma(alpha_vb[n,] + beta_vb[n,])
    E_log_1_h_n[,n] <- digamma(beta_vb[n,]) -  digamma(alpha_vb[n,] + beta_vb[n,])
  }
  
  # q_Z ===================================================
  for (n in 1:N){
    for(f in 1:F){
      E_Z[f,,n] <- exp(E_log_w[f,] + 
                         (V[f,n] * E_log_h[,n] + (1-V[f,n]) * E_log_1_h_n[,n]))
      E_Z[f,,n] <- E_Z[f,,n] / sum(E_Z[f,,n])
    }
  }

  # Compute lower bound ------------------------------------------------------
  if(TRUE){

      lower <- 0
      
      # E[log p(W | gamma)]
      lower <- lower + F*(lgamma(K*gamma) - K*lgamma(gamma)) - (gamma-1) * sum(E_log_w)
      
      # E[log p(Z | W)]
      lower <- lower + sum(sweep(E_Z, c(1,2), E_log_w, '*', check.margin=FALSE))
      
      # E[log p(H | alpha, beta)]      
      lower <- lower + 
               lgamma(alpha+beta) - lgamma(alpha) - lgamma(beta) +
               (alpha-1)*sum(E_log_h) + (beta-1)*sum(E_log_1_h_n)
      
      # E[log p(V | H, Z)]
      for(n in 1:N){
        for(f in 1:F){
          for(k in 1:K){
            lower <- lower + 
                     E_Z[f,k,n] * (V[f,n] * E_log_h[k,n] +
                                   (1-V[f,n]) * E_log_1_h_n[k,n])
          }
        }
      }
    
      # E[log q(W)]
      for(f in 1:F){
        lower <- lower + lgamma(sum(gamma_vb[f,]))
        for(k in 1:K){
          lower <- lower - lgamma(gamma_vb[f,k])
        }
        lower <- lower + sum(gamma_vb[f,]*E_log_w[f,])
      }
    
      # E[q(Z)]
      lower <- lower + sum(E_Z * log(E_Z))
      
      # E[q(H)]
      lower <- lower + sum(t(alpha_vb) * E_log_h + t(beta_vb) * E_log_1_h_n)
      
      
      #cat("\n Z colSumns:", apply(E_Z, 2, sum))
      cat("\n Lower bound:", lower)
  }
  # Store sample
  plot(sort(apply(E_Z, 2, sum)), main=j)
}


E_W <- compute_expectation_W_from_EZ(E_Z, gamma)
rownames(E_W) <- rownames(V)
plot_dictionary(E_W[,idx.sorted.W], Kmax=20, labels=TRUE)


