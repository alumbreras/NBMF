# This script is intended to familiriaze yourself with the
# general functions of the code.

# Generate a synthetic V matrix from a Dir-Ber model
# and then do two inferences assuming a Dir-Ber and a Dir-Dir model,
# respectively

devtools::load_all()
library(gtools)

# Generate synthetic data ----
alpha <- 0.1
beta <- 0.1
eta <- 0.1

F <- 100
N <- 100
K <- 4
W <- array(NA, dim = c(F,K))
H <- array(NA, dim = c(K,N))

for(f in 1:F){
  W[f,] <- rdirichlet(1, rep(alpha,K))
}
for(n in 1:N){
  H[,n] <- rbeta(K, alpha,beta)
}
V <- sample_V(W,H)

# Infer latent factors and posterior mean ----
K <-  100

# Assume a Dir-Ber model
res <- BernoulliNMF(V, K=K, model="DirBer", 
                                 alpha=1, beta=1, gamma=1/K, 
                                 iter=100, burnin = 0.9)
E_Z <- res$E_Z
E_W <- res$E_W
E_H <- res$E_H
E_V <- expectation(res)

plot_trace_size_components(res)
plot_dictionary(E_W, Kmax=10, rowlabels=TRUE, aspect.ratio = 7)
plot_V(E_V)
loglikelihood(res, V)$loglikelihood

# Assume a Dir-Dir model
res <- BernoulliNMF(V, K=K, model="DirDir", 
                                 alpha=1, gamma=1/K, 
                                 iter=100, burnin = 0.9)
E_Z <- res$E_Z
E_W <- res$E_W
E_H <- res$E_H
E_V <- expectation(res)

plot_trace_size_components(res)
plot_dictionary(E_W, Kmax=10, rowlabels=TRUE, aspect.ratio = 7)
plot_V(E_V)
loglikelihood(res, V)$loglikelihood
