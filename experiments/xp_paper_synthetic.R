# Create synthetic data from the three generative models:
# Beta-Dir
# Dir-Beta
# Dir-Dir 
# and learn their dictionaries

devtools::load_all()
library(gtools)

#' @title Plot matrix after applying hclust
#' @description Plot a matrix V after applying re-sorting its rows and columns
#' with a hclust algorithm in order to emphasize its structure
plot_V_hclust <- function(V){
  res <- heatmap(V)
  idx.row <- res$rowInd
  idx.col <- res$colInd
  V <- V[idx.row, idx.col]
  p <- plot_V(V)
  return(p)
}

F <- 100
N <- 100
K <- 4
W <- array(NA, dim = c(F,K))
H <- array(NA, dim = c(K,N))

# Generate non-structured matrices ----
alpha <- 1
beta <- 1
eta <- 1

# Beta-Dir
for(f in 1:F){
  W[f,] <- rdirichlet(1, rep(alpha,K))
}
for(n in 1:N){
  H[,n] <- rbeta(K, alpha,beta)
}
V <- sample_V(W,H)

plot_V(V)
p <- plot_V_hclust(V)
ggsave(p, 
       filename = paste0("fig_2_a_synthetic_BetaDir_nostruct.eps"), 
       height=5, 
       width=5, 
       units='cm')

# Dir-Beta
for(f in 1:F){
  W[f,] <- rbeta(K, alpha,beta)
}
for(n in 1:N){
  H[,n] <- rdirichlet(1, rep(alpha,K))
}

V <- sample_V(W,H)
plot_V(V)
p <- plot_V_hclust(V)
ggsave(p, 
       filename = paste0("fig_2_b_synthetic_DirBeta_nostruct.eps"), 
       height=5, 
       width=5, units='cm')

# Dir-Dir
for(f in 1:F){
  W[f,] <- rdirichlet(1, rep(alpha,K))
}
for(n in 1:N){
  H[,n] <- rdirichlet(1, rep(eta,K))
}
V <- sample_V(W,H)

plot_V(V)
p <- plot_V_hclust(V)
ggsave(p, 
       filename = paste0("fig_2_c_synthetic_DirDir_nostruct.eps"), 
       height=5, 
       width=5, 
       units='cm')

# Generate structured matrices ----
alpha <- 0.1
beta <- 0.1
eta <- 0.1

# Beta-Dir
for(f in 1:F){
  W[f,] <- rdirichlet(1, rep(alpha,K))
}
for(n in 1:N){
  H[,n] <- rbeta(K, alpha,beta)
}
V <- sample_V(W,H)

plot_V(V)
p <- plot_V_hclust(V)
ggsave(p, 
       filename = paste0("fig_2_a_synthetic_BetaDir_struct.eps"), 
       height=5, 
       width=5, 
       units='cm')

# Dir-Beta
for(f in 1:F){
  W[f,] <- rbeta(K, alpha,beta)
}
for(n in 1:N){
  H[,n] <- rdirichlet(1, rep(alpha,K))
}

V <- sample_V(W,H)
plot_V(V)
p <- plot_V_hclust(V)
ggsave(p, 
       filename = paste0("fig_2_b_synthetic_DirBeta_struct.eps"), 
       height=5, 
       width=5, units='cm')

# Dir-Dir
for(f in 1:F){
  W[f,] <- rdirichlet(1, rep(alpha,K))
}
for(n in 1:N){
  H[,n] <- rdirichlet(1, rep(eta,K))
}
V <- sample_V(W,H)

plot_V(V)
p <- plot_V_hclust(V)
ggsave(p, 
       filename = paste0("fig_2_c_synthetic_DirDir_struct.eps"), 
       height=5, 
       width=5, 
       units='cm')
