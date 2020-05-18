#
# Test sensitivity of concentration parameter over a toy dataset
# TODO: test for different datasizes. IN theory, with infinite data,
# we would reach zero sensitivity?

library(logisticPCA)
library(digest)
library(foreach)
library(doParallel)
library(bigmemory)
library(Matrix)
library(svs) # Classic NMFs
#library(NMF)
library(gtools)

devtools::load_all()

results_file = "results_sensitivity_criteo.csv"


gibbs.samples <- 5000
burnin <- 0.9

vb.iters      <- 500
repetitions   <- 3

df.results <- data.frame()
xp <-1

# A triple loop: repetitions, size, alpha
for(xp in 1:repetitions){
  
  for(size in seq(100, 400, by=100)){
  #for(size in c(1100, 1200)){
    F <- size
    N <- size
  
    K <- 10
    W <- array(NA, dim = c(F,K))
    H <- array(NA, dim = c(K,N))
    V.BetaDir <- array(NA, dim = c(F,N))
    
    # Beta-Dir
    alpha <- 1
    beta <- 1
    eta <- 1
    for(f in 1:F){
      W[f,] <- rdirichlet(1, rep(eta,K))
    }
    
    for(n in 1:N){
      H[,n] <- rbeta(K, alpha,beta)
    }
    
    
    for(f in 1:F){
      for(n in 1:N){
        V.BetaDir[f,n] <- rbinom(1, 1, prob = W[f,]%*%H[,n])
      }
    }
    
    V <- V.BetaDir
    
    
    for(gamma_var in seq(0.1, 10, by=1)){
      
      F <- dim(V)[1]
      N <- dim(V)[2]
    
      # Data -----------------------------------------------------------
      # Create a random training/test split
      # Insert some NA (these elements will be used for testing predictions)
      mask_test <- array(0, dim=c(F,N))
      for(f in 1:F){
        for(n in 1:N){
          if(runif(1) > 0.75){
            mask_test[f,n] = 1
          }
        }
      }
      # Test matrix
      V.test  <- as.matrix(V)
      is.na(V.test) <- !(as.logical(mask_test))
      
      # Training matrix
      V.train <- as.matrix(V)
      is.na(V.train) <- as.logical(mask_test)
      
      hash <- digest(V.train)
      ntest <- sum(mask_test)
      
      # Train -----------------------------------------------------------
      beta = 1
      gamma = 1
      #TODO K<10 throws index out of bounds exceptions
      modelGibbsDirBer <- BernoulliNMF(V.train, K=K, model="DirBer",
                                       alpha=alpha, beta=beta, gamma=gamma_var/K,
                                       iter=gibbs.samples, burnin = burnin)
      
      # Test -----------------------------------------------------------
      pred <- loglikelihood(modelGibbsDirBer, V.test)
      df.results <- list(xphash = hash,
                         xp=xp,
                         model= "Gibbs Beta-Dir",
                         loglikelihood = pred$loglikelihood,
                         K=K,
                         alpha=gamma_var,
                         ntest = ntest,
                         size=size)
      save_result(file = results_file, df.results)
    }
  }
}


#
# Plots ------------------------------------------------------------------------
#

# ******************************************************************************
# Prepare dataframe with properly named and sorted models
# ******************************************************************************

# Read all saved results from files
df <- read.table(file = results_file, sep = ",", header = T) %>% distinct()

# compute total perplexity
df <- df %>% mutate(perplexity = -loglikelihood/ntest)

# ******************************************************************************
# Plot from dataframes
# ******************************************************************************

# Perplexity (used in the paper)
base_size <- 8
p <- ggplot(df, aes(x=alpha, y=perplexity, group=alpha)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(colour = "red", size = 1, shape=3) +
  facet_grid(size ~ ., scales = "fixed") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, size=base_size*1, hjust = 1, colour = "black"),
        axis.text.y = element_text(size= base_size, color = 'black'),
        strip.background =element_rect(fill="white"),
        aspect.ratio = 1/3)
print(p)
ggsave(p, filename = "fig_boxplots_sensitivity_predictions_perplexity2.eps", height=18, width=13, units='cm')
