# Easy, synthetic example with a block diagonal matrix to test the 
# self-regularization property of the model
# Gibbs sampling with p(z_i | Z_{-i}) and marginalizing W

library(dplyr)
library(tidyr)
library(ggplot2)
library(Matrix)
library(gtools)


x <- list(matrix(1,ncol=6, nrow=6), 
          matrix(1,ncol=6, nrow=6),
          matrix(1,ncol=6, nrow=6))
V <- as.matrix(bdiag(x))



# Hypersimple
V <- diag(1,3)

# High F
x <- list(matrix(1,nrow=50, ncol=6), 
          matrix(1,nrow=50, ncol=6),
          matrix(1,nrow=50, ncol=6))
V <- as.matrix(bdiag(x))


x <- list(matrix(1,nrow=6, ncol=50), 
          matrix(1,nrow=6, ncol=50),
          matrix(1,nrow=6, ncol=50))
V <- as.matrix(bdiag(x))


x <- list(matrix(1,nrow=50, ncol=50), 
          matrix(1,nrow=50, ncol=50),
          matrix(1,nrow=50, ncol=50))
V <- as.matrix(bdiag(x))


load("../data/adjacency_parlamentcat_sorted.RData")

nsamples <- 100

F <- dim(V)[1]
N <- dim(V)[2]
K <- 40
alpha <- 1/K
beta <- 1/K
gamma <- 0.01
gamma <- gamma/K
initW <- matrix(rbeta(F*K, shape1=1, shape2=1), F,K)
W <- initW
#initW <- W*0
#initW[,1] <- 1
rownames(W) <- rownames(V)


# Gibbs ------------------------------------------------------------------------
F <- dim(W)[1]
K <- dim(W)[2]

Z_samples <- array(NA, dim=c(nsamples, F, K, N))
Z <- array(0, dim=c(F, K, N))
for(n in 1:N){
  for(f in 1:F){
    Z[f,,n] <- c(rmultinom(1, 1, prob=rep(1,K)))
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
    v_n <- V[,n]
    choices <- rep(NA, F)
    for(f in 1:F){
      Z[f,,n] <- 0 # remove current value
      L_plus <- colSums(Z[,,n]*v_n) # sum over f
      L_minus <- colSums(Z[,,n]*(1-v_n)) # sum over f
      Z_minus_n <- rowSums(Z[f,,])
      probs0 <- v_n[f] * (alpha + L_plus) + (1-v_n[f]) * (beta + L_minus)
      probs1 <- probs0 / (alpha + beta + L_plus + L_minus)
      probs2 <- (gamma + Z_minus_n) 
      probs <- (probs1 * probs2)
      probs <- probs /sum(probs)
      k_chosen <- sample(K, size=1, prob=probs) # assign new component to f
      Z[f,k_chosen,n] <- 1 # re-assign value to chosen k
      choices[f] <- k_chosen
    }
  }
  cat("\n Z colSumns:", colSums(rowSums(Z, dims = 2)))
  # Store sample
  Z_samples[j,,,] <- Z
  plot(sort(colSums(rowSums(Z, dims = 2))))
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

