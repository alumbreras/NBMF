# Easy, synthetic example with a block diagonal matrix to test the 
# self-regularization property of the model
# Gibbs sampling with p(z_i | Z_{-i}) and marginalizing W

library(dplyr)
library(tidyr)
library(ggplot2)
library(Matrix)
library(gtools)

load("../data/adjacency_parlamentcat_sorted.RData")

nsamples <- 1000

F <- dim(V)[1]
N <- dim(V)[2]
K <- 50
alpha <- 1
gamma <- 1/K


# Gibbs ------------------------------------------------------------------------

Z_samples <- array(NA, dim=c(nsamples, F, K, N))
C_samples <- array(NA, dim=c(nsamples, F, K, N))
Z <- array(0, dim=c(F, K, N))
C <- array(0, dim=c(F, K, N))

for(n in 1:N){
  for(f in 1:F){
    Z[f,,n] <- c(rmultinom(1, 1, prob=rep(1,K)))
    C[f,,n] <- c(rmultinom(1, 1, prob=rep(1,K)))
  }
}
#Z[,1,] <- V # all asssignments to cluster 1

j <- 1
n <- 1
f <- 1
for(j in 1:nsamples){
  if(!j %% 20){cat('\n sample:', j)}
  cat("\n*******************************")
  for (n in 1:N){
    for(f in 1:F){
      
      # remove current values and compute sums
      Z[f,,n] <- 0 
      C[f,,n] <- 0
      C_sum <- colSums(C[,,n])
      Z_sum <- rowSums(Z[f,,])
      
      if(V[f,n]){
        probs <- (gamma + Z_sum)*(alpha + C_sum)
        probs <- probs / sum(probs)
        k_chosen <- sample(K, size=1, prob=probs) # assign new component to f
        Z[f,k_chosen,n] <- 1 # re-assign value to chosen k
        C[f,k_chosen,n] <- 1 # re-assign value to chosen k
      }else{
        probs_z <- (gamma + Z_sum)*(1-C[f,,n]) # BUG:: C[f,,n] has ben removed
        probs_z <- probs_z / sum(probs_z)
        k_chosen_z <- sample(K, size=1, prob=probs_z) # assign new component to f
        Z[f,k_chosen_z,n] <- 1 # re-assign value to chosen k
        
        probs_c <- (alpha + C_sum)*(1-Z[f,,n])
        probs_c <- probs_c / sum(probs_c)
        k_chosen_c <- sample(K, size=1, prob=probs_c) # assign new component to f
        C[f,k_chosen_c,n] <- 1 # re-assign value to chosen k
      }
    }
  }

  # Store samples
  Z_samples[j,,,] <- Z
  C_samples[j,,,] <- C
  
  #cat("\n Z colSumns:", colSums(rowSums(Z, dims = 2)))
  #cat("\n C colSumns:", colSums(rowSums(C, dims = 2)))
  plot(sort(colSums(rowSums(Z, dims = 2))))
  #plot(sort(colSums(rowSums(C, dims = 2))))  
}

# Estimate W -------------------------------------------------------------------
E_Z <- colMeans(Z_samples[1:j-1,,,])      # average over the J samples
E_Z_sum <- rowSums(E_Z, dims=2)
E_W <- array(NA, dim = c(F,K))
for(f in 1:F){
  E_W[f,] = gamma + E_Z_sum[f,]
}

for(f in 1:F){
  E_W[f,] <- E_W[f,] /sum(E_W[f,])
}
W <- E_W

# Plot -------------------------------------------------------------------------
W <- W[,order(-colSums(W))]
rownames(W) <- 1:dim(W)[1]
colnames(W) <- 1:dim(W)[2]

# Plot final dictionary W ------------------------------------------------------
df.W <- data.frame(W)
colnames(df.W) <- 1:dim(W)[2]
df.W$name <- 1:dim(W)[1]

df.W <- gather(df.W, key='dimension', value='value', -name)
df.W$dimension <- as.numeric(df.W$dimension)

text_components <-  10
text_features <- 6
p <- ggplot(df.W, aes(x=factor(dimension), y=name)) + 
  geom_tile(aes(fill = value), colour = "white") + 
  scale_fill_gradient(low = "white", high = "steelblue", na.value = 'red') +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  xlab("component") +
  theme_bw() +
  theme(axis.text.x = element_text(size=text_components,  colour = "grey50"),
        axis.text.y = element_text(size=text_features, hjust = 1, colour = "grey50"),
        axis.ticks.y = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")
print(p)

