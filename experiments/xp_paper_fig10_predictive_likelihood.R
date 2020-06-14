devtools::load_all()
library(logisticPCA)
library(digest)
library(foreach)
library(doParallel)
library(bigmemory)
library(Matrix)
library(svs) # Classic NMFs

results_file = "results_predictions_test.csv"
plots_file <- "fig10_boxplots_predictions_perplexity.eps"

# Too big. Just show the dictionary
#load("../data/unvotes100.RData")
data("animals")
data("paleo")
data("catalanparliament")
data("unvotes100coldwar_absna")
unvotes100coldwar <- unvotes100coldwar[,1:2500]

# Datasets used in this experiment in the paper
dataset_names <- c('animals', 'lastfm', 'paleo', 'parlamentcat', 'unvotes')
datasets <- list(animals, lastfm, paleo, catalanparliament, unvotes100coldwar)

# Datasets to do a fast check
dataset_names <- c("animals", "parlamentcat")
datasets <- list(animals, catalanparliament)

# Training parameters
# Use small parameters for fast checks 
# Use larger parameters for accurate results 
gibbs.samples <- 5
burnin <- 0.9
vb.iters <- 10
repetitions   <- 3

models.BetaDirGS <- TRUE
models.BetaDirVB <- TRUE
models.DirBerVB  <- TRUE
models.DirDirGS  <- TRUE
models.cbICA     <- TRUE
models.bICA      <- TRUE
models.logPCA    <- TRUE

for (i in 1:length(dataset_names)) {
  dataset <- dataset_names[i]
  V <- datasets[[i]]
  F <- dim(V)[1]
  N <- dim(V)[2]
  
  for(xp in 1:repetitions){
    df.results <- data.frame()
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
    V.test  <- as.matrix(V)
    V.train <- as.matrix(V)
  
    # Test matrix
    is.na(V.test) <- !(as.logical(mask_test))
    # Training matrix
    is.na(V.train) <- as.logical(mask_test)
  
    hash <- digest(V.train)
    ntest <- sum(mask_test)
  
    # Gibbs DirBer --------------------------------------------------------------------
    if(models.BetaDirGS){
      K = 100
      alpha = 1
      beta = 1
      gamma = 1
      modelGibbsDirBer <- BernoulliNMF(V.train, 
                                       K = K, 
                                       model = "DirBer",
                                       alpha = alpha, 
                                       beta = beta, 
                                       gamma = gamma/K,
                                       iter = gibbs.samples, 
                                       burnin = burnin)
      E_Z <- modelGibbsDirBer$E_Z
      E_W <- modelGibbsDirBer$E_W
      E_H <- modelGibbsDirBer$E_H
      pred <- loglikelihood(modelGibbsDirBer, V.test)
      df.results <- list(date = date(),
                         xphash = hash,
                         xp = xp,
                         dataset = dataset,
                         ntest = ntest,
                         model= "Gibbs Beta-Dir",
                         K = K,
                         loglikelihood = pred$loglikelihood
                         )
      save_result(file = results_file, df.results)
    }
    
    # Gibbs DirDir ---------------------------------------------------------------
    if(models.DirDirGS){
      K = 100
      alpha = 1
      beta = 1
      gamma = 1
      modelGibbsDirDir <- BernoulliNMF(V.train, 
                                       K = K, 
                                       model ="DirDir",
                                       alpha = alpha, 
                                       gamma = gamma/K,
                                       iter = gibbs.samples, 
                                       burnin = burnin)
      E_Z <- modelGibbsDirDir$E_Z
      E_W <- modelGibbsDirDir$E_W
      E_H <- modelGibbsDirDir$E_H
      pred <- loglikelihood(modelGibbsDirDir, V.test)
      df.results <- list(date = date(),
                         xphash = hash,
                         xp = xp,
                         dataset = dataset,
                         ntest = ntest,
                         model = "Gibbs Dir-Dir",
                         K = K,
                         loglikelihood = pred$loglikelihood
                         )
      save_result(file = results_file, df.results)
    }
    
    # VB DirBer ------------------------------------------------------------------
    if(models.DirBerVB){
      K = 100
      alpha = 1
      beta = 1
      gamma = 1
      modelVBDirBer <- BernoulliNMF(V.train, 
                                    K = K, 
                                    model = "DirBer", 
                                    method = "VB",
                                    alpha = alpha, 
                                    beta = beta, 
                                    gamma = gamma/K,
                                    iter = vb.iters)
      
      E_Z <- modelVBDirBer$E_Z
      E_W <- modelVBDirBer$E_W
      E_H <- modelVBDirBer$E_H
      pred <- loglikelihood(modelVBDirBer, V.test)
      df.results <- list(date = date(),
                         xphash = hash,
                         xp = xp,
                         dataset = dataset,
                         ntest = ntest,
                         model = "VB Beta-Dir",
                         K = K,
                         loglikelihood = pred$loglikelihood)
      save_result(file = results_file, df.results)
    }
    
    for(k in 2:7){
      # VB Aspect --------------------------------------------------------------
      if(models.bICA){
      alpha = 1
      beta = 1
      gamma = 1
      modelVBaspect <- BernoulliNMF(V.train, 
                                    K=k, 
                                    model="aspectRcpp", 
                                    method="VB",
                                    alpha=alpha, 
                                    beta=beta, 
                                    gamma=gamma,
                                    iter=vb.iters)
      E_Z <- modelVBaspect$E_Z
      E_W <- modelVBaspect$E_W
      E_H <- modelVBaspect$E_H
      pred <- loglikelihood(modelVBaspect, V.test)
      df.results <- list(date = date(),
                         xphash = hash,
                         xp = xp,
                         dataset = dataset,
                         ntest = ntest,
                         model = "Aspect",
                         K = k,
                         loglikelihood = pred$loglikelihood)
  
      }
      # collapsed Aspect-------------------------------------------------------------
      if(models.cbICA){
        modelVBaspectcol <- BernoulliNMF(V.train, K=k, model="DirBer", method="VB",
                                         alpha=alpha, beta=beta, gamma=gamma,
                                         iter=vb.iters)
        E_Z <- modelVBaspectcol$E_Z
        E_W <- modelVBaspectcol$E_W
        E_H <- modelVBaspectcol$E_H
        pred <- loglikelihood(modelVBaspectcol, V.test)
        df.results <- list(date = date(),
                           xphash = hash,
                           xp = xp,
                           dataset = dataset,
                           ntest = ntest,
                           model = "cAspect",
                           K = k,
                           loglikelihood = pred$loglikelihood)
        save_result(file = results_file, df.results)
      }
      
      # Logistic PCA-------------------------------------------------------------
      if(models.logPCA){
        logpca_model <- logisticPCA::logisticSVD(V.train, k, max_iters = 2000)
        W <- logpca_model$A
        H <- logpca_model$B
        E_V <- fitted(logpca_model, type = "response")
        V.probs <- V.test*E_V + (1-V.test)*(1-E_V)
        like <- sum(log(V.probs), na.rm = TRUE)
        df.results <- list(date = date(),
                           xphash = hash,
                           xp = xp,
                           dataset = dataset,
                           ntest = ntest,
                           model = "logPCA",
                           K = k,
                           loglikelihood = like)
        save_result(file = results_file, df.results)
      }
    }
  }
}

# Plots ------------------------------------------------------------------------

# ******************************************************************************
# Prepare dataframe with properly named and sorted models
# ******************************************************************************

# Read all saved results from files
df <- read.table(file = results_file, sep = ",", header = T) %>% distinct()

# filter logPCA with high k (uninteresting)
df <- df[!(df$model == "logPCA" & df$K > 4),]
df <- df %>% filter(loglikelihood != -Inf)

# Pretty names. Denote each Aspect/bICA model as
# bICA-k where k in the number of components

df$model <- as.character(df$model)
df$model[df$model == "VB Beta-Dir"]    <- "1.*Beta-Dir VB"
df$model[df$model == "Gibbs Beta-Dir"] <- "2.*Beta-Dir GS"
df$model[df$model == "Gibbs Dir-Dir"]  <- "3.*Dir-Dir GS"
df$model[df$model == "cAspect"]        <- "4.*c-bICA"
df$model[df$model == "Aspect"]         <- "5.bICA"
df$model[df$model == "logPCA"]         <- "6.logPCA"

df$dataset <- as.character(df$dataset)
df$dataset[df$dataset == "unvotes100"]   <- "unvotes"
df$dataset[df$dataset == "parlament"]    <- "parliament"
df$dataset <- as.factor(df$dataset)

for(i in 1:nrow(df)){
  if(df$model[i] == "4.*c-bICA") {
    df$model[i] <- paste0("4.*c-bICA-", df$K[i])
  }
}
for(i in 1:nrow(df)){
  if(df$model[i] == "5.bICA") {
    df$model[i] <- paste0("5.bICA-", df$K[i])
  }
}
for(i in 1:nrow(df)){
  if(df$model[i] == "6.logPCA") {
    df$model[i] <- paste0("6.logPCA-", df$K[i])
  }
}

levels <- substring(sort(unique(df$model)), 3)
df$model <- substring(df$model,3)
df$model <- factor(df$model, levels = levels)

# compute total perplexity and
# loglikelihood with respect to the best score in each dataset
df <- df %>% mutate(perplexity = -loglikelihood/ntest)
df <- df %>% group_by(dataset) %>%
  mutate(rel_loglikelihood = loglikelihood - max(loglikelihood))

# ******************************************************************************
# Plot from dataframes
# ******************************************************************************

# Perplexity (used in the paper)
base_size <- 8
p <- ggplot(df,
            aes(x=model, y=perplexity)) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(dataset ~ ., scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, 
                                   size = base_size, 
                                   hjust = 1, 
                                   colour = "black"),
        axis.text.y = element_text(size= base_size, color = 'black'),
        strip.background = element_rect(fill = "white"),
        aspect.ratio = 1/4)
print(p)
ggsave(p, filename = plots_file, height = 18, width = 13, units = 'cm')