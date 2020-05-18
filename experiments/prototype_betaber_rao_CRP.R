#TODO: use only K arrays, adding and removing accordingly
load("../data/adjacency_parlamentcat_sorted.RData")

nsamples <- 100
F <- dim(V)[1]
N <- dim(V)[2]
K <- 20 # maximum number of components
alpha <- 1
beta <- 1
gamma <- 1
gamma <- gamma/K

Z_samples <- array(NA, dim=c(nsamples, F, K, N))
Z <- array(0, dim=c(F, K, N))
for(n in 1:N){
  for(f in 1:F){
    Z[f,,n] <- c(rmultinom(1, 1, prob=rep(1,K)))
  }
}


j <- 1
n <- 1
f <- 1
for(j in 1:nsamples){
  if(!j %% 5){cat('\n sample:', j)}
  cat("\n*******************************")
  for (n in 1:N){
    #cat("\nn:", n)
    v_n <- V[,n]
    n_pos <- colSums(Z[,,n]*v_n)
    n_neg <- colSums(Z[,,n]*(1-v_n))
    
    choices <- rep(NA, F)
    for(f in 1:F){
      Z[f,,n] <- 0 # remove current value
      
      # Utils
      L_plus <- colSums(Z[,,n]*v_n) # sum over f
      L_minus <- colSums(Z[,,n]*(1-v_n)) # sum over f
      Z_minus_n <- rowSums(Z[f,,]) # TODO: bottleneck
            
      # Prior and likelihood
      probs_CRP <- Z_minus_n
      probs_likelihood <- v_n[f] * (alpha + L_plus) + (1-v_n[f]) * (beta + L_minus)
      probs_likelihood <- probs_likelihood / (alpha + beta + L_plus + L_minus)
      components <- which(probs_CRP>0)
      first_empty_component <- which(Z_minus_n==0)[1]
      probs <- c(probs_CRP * probs_likelihood)
      probs <- probs[components] # a compact vector without holes
      probs <- c(probs, gamma/2) # add probability for empty component
      probs <- probs /sum(probs)                          # normalize
      idx <- c(components, first_empty_component)
      
      # Sample 
      k_chosen <- sample(idx, size=1, prob=probs)           
      Z[f,k_chosen,n] <- 1
    }
  }
  K <- sum(colSums(rowSums(Z, dims = 2)) > 0)
  cat('\nK = ', K)
  # Store sample
  Z_samples[j,,,] <- Z
  plot(sort(colSums(rowSums(Z, dims = 2))))
}


# Estimate W
E_Z <- colMeans(Z_samples[1:j-1,,,])      # average over the J samples
E_Z_sum <- rowSums(E_Z, dims=2)
E_W <- array(NA, dim = c(F,K))
E_W <- array(NA, dim = c(F,K+1))
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


