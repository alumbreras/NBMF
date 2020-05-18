# Easy, synthetic example with a block diagonal matrix to test the 
# self-regularization property of the model

library(dplyr)
library(tidyr)
library(ggplot2)
library(Matrix)
library(gtools)

x <- list(matrix(1,nrow=50, ncol=6), 
          matrix(1,nrow=50, ncol=6),
          matrix(1,nrow=50, ncol=6))
V <- as.matrix(bdiag(x))

x <- list(matrix(1,ncol=6, nrow=6), 
          matrix(1,ncol=6, nrow=6),
          matrix(1,ncol=6, nrow=6))
V <- as.matrix(bdiag(x))

x <- list(matrix(1,nrow=6, ncol=50), 
          matrix(1,nrow=6, ncol=50),
          matrix(1,nrow=6, ncol=50))
V <- as.matrix(bdiag(x))

#load("../data/adjacency_parlamentcat.RData")


nsamples <- 20

F <- dim(V)[1]
N <- dim(V)[2]
K <- 3
alpha <- 1
beta <- 1
alpha <- alpha/K
beta <- beta/K
alpha <- 1
beta <- 1
gamma <- 1
initW <- matrix(rbeta(F*K, shape1=1, shape2=1), F,K)
W <- initW
initW <- W*0
initW[,1] <- 1
rownames(W) <- rownames(V)


# Gibbs ------------------------------------------------------------------------
F <- dim(W)[1]
K <- dim(W)[2]

W_samples <- array(NA, dim=c(nsamples, F, K))
Z <- array(0, dim=c(F, K, N))
Z[,1,] <- V # all asssignments to cluster 1
W <- initW

j <- 1
n <- 1
f <- 1
for(j in 1:nsamples){
  if(!j %% 20){cat('\n sample:', j)}
  
  for (n in 1:N){
    v_n <- V[,n]
    for(f in 1:F){
      Z[f,,n] <- 0 # remove current value
      L_plus <- colSums(Z[,,n]*v_n) # sum over f
      L_minus <- colSums(Z[,,n]*(1-v_n)) # sum over f
      probs <- W[f,] * (alpha + L_plus) * (beta + L_minus) / (alpha + beta + L_plus + L_minus)
      probs <- probs/sum(probs)
      cat('\n probs', probs)
      k_chosen <- sample(K, size=1, prob=probs) # assign new component to f
      Z[f,k_chosen,n] <- 1 # re-assign value to chosen k
    }
  }
  
  # Sample W
  Z_sum <- rowSums(Z, dims=2)
  for(f in 1:F){
    gamma_post = gamma + Z_sum[f,]
    W[f,] <- rdirichlet(1, gamma_post)
  }
  #cat('\n W columns:', colSums(W))
  
  # Store sample
  W_samples[j,,] <- W
}
