library(logisticPCA)
library(digest)
library(foreach)
library(doParallel)
library(bigmemory)
library(Matrix)
library(svs) # Classic NMFs
library(NMF)

devtools::load_all()

results_file = "results_sensitivity_realdatasets.csv"

# Too big. Just show the dictionary
#load("../data/unvotes100.RData")
data("unvotes100coldwar")
V <- unvotes100coldwar[,1:2500]
dataset <- 'unvotes'

data("catalanparliament")
V <- catalanparliament
dataset <- "parlament"

data("unvotes100coldwar")
V <- unvotes100coldwar[,1:2500]
dataset <- 'unvotes'

data("catalanparliament")
V <- catalanparliament
dataset <- "parlament"

data("animals")
V <- animals
dataset <- "animals"

data("lastfm")
V <- lastfm
dataset <- "lastfm"

data(paleo)
V <- paleo
dataset <- "paleo"


models.BetaDirGS <- TRUE
models.BetaDirVB <- TRUE
models.DirDir    <- TRUE
models.cbICA     <- TRUE
models.bICA      <- TRUE
models.logPCA    <- TRUE
models.KL        <- TRUE
models.NMF       <- TRUE

F <- dim(V)[1]
N <- dim(V)[2]

gibbs.samples <- 5000
burnin <- 0.9

vb.iters      <- 500
repetitions   <- 3

for(alpha_var in c(0.1, 1, 2, 3, 4, 5, 6, 7, 8, 9)){
for(xp in 1:repetitions){
  df.results <- data.frame()
  
  F <- dim(V)[1]
  N <- dim(V)[2]
  
  #V.big <- attach.big.matrix(V.desc)
  
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
  #is.na(V.test) <- !as.matrix(as.logical(mask_test))
  is.na(V.test) <- !(as.logical(mask_test))
  
  # Training matrix
  #is.na(V.train) <- as.matrix(as.logical(mask_test))
  is.na(V.train) <- as.logical(mask_test)
  
  
  hash <- digest(V.train)
  ntest <- sum(mask_test)
  
  # Gibbs DirBer --------------------------------------------------------------------
  if(models.BetaDirGS){
    K = 100
    beta = 1
    gamma = 1
    modelGibbsDirBer <- BernoulliNMF(V.train, K=K, model="DirBer",
                                     alpha=alpha_var, beta=beta, gamma=gamma/K,
                                     iter=gibbs.samples, burnin = burnin)
    E_Z <- modelGibbsDirBer$E_Z
    E_W <- modelGibbsDirBer$E_W
    E_H <- modelGibbsDirBer$E_H
    pred <- loglikelihood(modelGibbsDirBer, V.test)
    df.results <- list(xphash = hash,
                       xp=xp,
                       model= "Gibbs Beta-Dir",
                       loglikelihood = pred$loglikelihood,
                       K=K,
                       alpha = alpha_var,
                       dataset = dataset,
                       ntest = ntest)
    save_result(file = results_file, df.results)
  }
  
  # Gibbs DirDir ---------------------------------------------------------------
  if(models.DirDirGS){
    K = 100
    beta = 1
    gamma = 1
    modelGibbsDirDir <- BernoulliNMF(V.train, K=K, model="DirDir",
                                     alpha=alpha_var, gamma=gamma/K,
                                     iter=gibbs.samples, burnin = burnin)
    E_Z <- modelGibbsDirDir$E_Z
    E_W <- modelGibbsDirDir$E_W
    E_H <- modelGibbsDirDir$E_H
    pred <- loglikelihood(modelGibbsDirDir, V.test)
    df.results <- list(xphash = hash,
                       xp = xp,
                       model= "Gibbs Dir-Dir",
                       loglikelihood = pred$loglikelihood,
                       K=K,
                       alpha = alpha_var,
                       dataset = dataset,
                       ntest = ntest)
    save_result(file = results_file, df.results)
  }
  
  # VB DirBer ------------------------------------------------------------------
  if(models.DirBerVB){
    K = 100
    alpha = 1
    beta = 1
    gamma = 1
    modelVBDirBer <- BernoulliNMF(V.train, K=K, model="DirBer", method="VB",
                                  alpha=alpha, beta=beta, gamma=gamma/K,
                                  iter=vb.iters)
    
    E_Z <- modelVBDirBer$E_Z
    E_W <- modelVBDirBer$E_W
    E_H <- modelVBDirBer$E_H
    pred <- loglikelihood(modelVBDirBer, V.test)
    df.results <- list(xphash = hash,
                       xp = xp,
                       model= "VB Beta-Dir",
                       loglikelihood = pred$loglikelihood,
                       K=K,
                       alpha = alpha_var,
                       dataset = dataset,
                       ntest = ntest)
    save_result(file = results_file, df.results)
  }
  
  # VB Aspect ----------------------------------------------------------------
  for(k in 2:7){
    if(models.bICA){
      beta = 1
      gamma = 1
      modelVBaspect <- BernoulliNMF(V.train, K=k, model="aspectRcpp", method="VB",
                                    alpha=alpha_var, beta=beta, gamma=gamma,
                                    iter=vb.iters)
      E_Z <- modelVBaspect$E_Z
      E_W <- modelVBaspect$E_W
      E_H <- modelVBaspect$E_H
      pred <- loglikelihood(modelVBaspect, V.test)
      df.results <- list(xphash = hash,
                        xp = xp,
                        model= "Aspect",
                        loglikelihood = pred$loglikelihood,
                        K=k,
                        alpha = alpha_var,
                        dataset = dataset,
                        ntest = ntest)
      save_result(file = results_file, df.results)
    }
    
    # collapsed Aspect-------------------------------------------------------------
    if(models.cbICA){
      modelVBaspectcol <- BernoulliNMF(V.train, K=k, model="DirBer", method="VB",
                                       alpha=alpha_var, beta=beta, gamma=gamma,
                                       iter=vb.iters)
      
      E_Z <- modelVBaspectcol$E_Z
      E_W <- modelVBaspectcol$E_W
      E_H <- modelVBaspectcol$E_H
      pred <- loglikelihood(modelVBaspectcol, V.test)
      df.results <- list(xphash = hash,
                         xp = xp,
                         model= "cAspect",
                         loglikelihood = pred$loglikelihood,
                         K=k,
                         alpha = alpha_var,
                         dataset = dataset,
                         ntest = ntest)
      save_result(file = results_file, df.results)
    }
    
    # NMF-KL-------------------------------------------------------------
    if(models.KL){
      # NMF Lee and Seung (KL/Poisson)
      #modelNMF_KL <- NMF_KL(V.train, K=k, alpha=alpha, beta=beta, gamma=gamma, 
      #                      iter=vb.iters)
      modelNMF_KL <- NMF::nmf(V, k, "KL")
      W <- basis(modelNMF_KL)
      H <- coef(modelNMF_KL)
      E_V <- W %*% H
      V.probs <- V.test*E_V + (1-V.test)*(1-E_V)
      like <- sum(log(V.probs), na.rm = TRUE)
      df.results <- list(xphash = hash,
                         xp = xp,
                         model= "KL",
                         loglikelihood = like,
                         K=k,
                         alpha = alpha_var,
                         dataset = dataset,
                         ntest = ntest)
      save_result(file = results_file, df.results)
    }
    
    # Logistic PCA -------------------------------------------------------------
    if(models.logPCA){
      logpca_model <- logisticPCA::logisticSVD(V.train, k, max_iters = 2000)
      W <- logpca_model$A
      H <- logpca_model$B
      E_V <- fitted(logpca_model, type = "response")
      V.probs <- V.test*E_V + (1-V.test)*(1-E_V)
      like <- sum(log(V.probs), na.rm = TRUE)
      df.results <- list(xphash = hash,
                         xp = xp,
                         model= "logPCA",
                         loglikelihood = like,
                         K=k,
                         alpha = alpha_var,
                         dataset = dataset,
                         ntest = ntest)
      save_result(file = results_file, df.results)
    }
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

# filter logPCA with high k (uninteresting)
df <- df[!(df$model == "logPCA" & df$k > 4),]
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
df$model[df$model == "KL"]             <- "7.KL"

df$dataset <- as.character(df$dataset)
df$dataset[df$dataset == "unvotes100"]   <- "unvotes"
df$dataset[df$dataset == "parlament"]    <- "parliament"
df$dataset <- as.factor(df$dataset)

for(i in 1:nrow(df)){
  if(df$model[i] == "4.*c-bICA") {
    df$model[i] <- paste0("4.*c-bICA-", df$k[i])
  }
}
for(i in 1:nrow(df)){
  if(df$model[i] == "5.bICA") {
    df$model[i] <- paste0("5.bICA-", df$k[i])
  }
}
for(i in 1:nrow(df)){
  if(df$model[i] == "6.logPCA") {
    df$model[i] <- paste0("6.logPCA-", df$k[i])
  }
}
for(i in 1:nrow(df)){
  if(df$model[i] == "7.KL") {
    df$model[i] <- paste0("7.KL-", df$k[i])
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

# Loglikelihood with whiskers and crosses
base_size <- 8
p <- ggplot(df %>% filter(!grepl("KL", model)), aes(x=model, y=rel_loglikelihood)) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(dataset ~ ., scales = "free") +
  geom_point(colour = "red", size = 1, shape=3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, size=base_size*1, hjust = 1, colour = "black"),
        aspect.ratio = 1/5) + ylab("loglikelihood")
print(p)
ggsave(p, filename = "fig_boxplots_predictions_loglikelihood.eps",
       height=8, width=8, units='cm')

# Perplexity (used in the paper)
base_size <- 8
p <- ggplot(df %>% filter(!grepl("KL", model), dataset != "epinions100"),
            aes(x=model, y=perplexity)) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(dataset ~ ., scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, size=base_size*1, hjust = 1, colour = "black"),
        axis.text.y = element_text(size= base_size, color = 'black'),
        strip.background =element_rect(fill="white"),
        aspect.ratio = 1/4)
print(p)
ggsave(p, filename = "fig_boxplots_predictions_perplexity.eps", height=18, width=13, units='cm')
