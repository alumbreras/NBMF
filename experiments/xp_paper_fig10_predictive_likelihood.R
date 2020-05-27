library(logisticPCA)
library(digest)
library(foreach)
library(doParallel)
library(bigmemory)
library(Matrix)
library(svs) # Classic NMFs
#library(NMF)
library(lda) # LDA
library(topicmodels) # LDA (TODO install it and use perplexity function)

# TODO: lastfm train may have all zeros and that be a problem for LDA?
# What is the effect of this for the other methods?

results_file = "results_predictions_criteo.csv"
devtools::load_all()


data("catalanparliament")
V <- catalanparliament
dataset <- "parlament"

data("paleo")
V <- paleo
dataset <- "paleo"

data("catalanparliament")
V <- catalanparliament
dataset <- "parlament"

# Too big. Just show the dictionary
#load("../data/unvotes100.RData")
data("unvotes100coldwar")
V <- unvotes100coldwar[,1:2500]
dataset <- 'unvotes'

data("lastfm")
V <- lastfm
dataset <- "lastfm"

data("animals")
V <- animals
dataset <- "animals"


F <- dim(V)[1]
N <- dim(V)[2]

gibbs.samples <- 5000
burnin <- 0.9

vb.iters      <- 500
repetitions   <- 5

models.BetaDirGS <- FALSE
models.BetaDirVB <- FALSE
models.DirBerVB  <- FALSE
models.DirDirGS  <- FALSE
models.cbICA     <- FALSE
models.bICA      <- FALSE
models.logPCA    <- FALSE
models.KL        <- FALSE
models.NMF       <- FALSE
models.LDA       <- TRUE

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
  #plot_V(mask_test)

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
    modelGibbsDirBer <- BernoulliNMF(V.train, K=K, model="DirBer",
                                     alpha=alpha, beta=beta, gamma=gamma/K,
                                     iter=gibbs.samples, burnin = burnin)
    E_Z <- modelGibbsDirBer$E_Z
    E_W <- modelGibbsDirBer$E_W
    E_H <- modelGibbsDirBer$E_H
    pred <- loglikelihood(modelGibbsDirBer, V.test)
    df.results <- list(date = date(),
                       xphash = hash,
                       xp=xp,
                       dataset = dataset,
                       ntest = ntest,
                       model= "Gibbs Beta-Dir",
                       K=K,
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
    modelGibbsDirDir <- BernoulliNMF(V.train, K=K, model="DirDir",
                                     alpha=alpha, gamma=gamma/K,
                                     iter=gibbs.samples, burnin = burnin)
    E_Z <- modelGibbsDirDir$E_Z
    E_W <- modelGibbsDirDir$E_W
    E_H <- modelGibbsDirDir$E_H
    pred <- loglikelihood(modelGibbsDirDir, V.test)
    df.results <- list(date = date(),
                       xphash = hash,
                       xp=xp,
                       dataset = dataset,
                       ntest = ntest,
                       model= "Gibbs Dir-Dir",
                       K=K,
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
    modelVBDirBer <- BernoulliNMF(V.train, K=K, model="DirBer", method="VB",
                                  alpha=alpha, beta=beta, gamma=gamma/K,
                                  iter=vb.iters)
    
    E_Z <- modelVBDirBer$E_Z
    E_W <- modelVBDirBer$E_W
    E_H <- modelVBDirBer$E_H
    pred <- loglikelihood(modelVBDirBer, V.test)
    df.results <- list(date = date(),
                       xphash = hash,
                       xp=xp,
                       dataset = dataset,
                       ntest = ntest,
                       model= "VB Beta-Dir",
                       loglikelihood = pred$loglikelihood)
    save_result(file = results_file, df.results)
  }
  
  for(k in 2:7){
    
    # VB Aspect ----------------------------------------------------------------
    if(models.bICA){
    alpha = 1
    beta = 1
    gamma = 1
    modelVBaspect <- BernoulliNMF(V.train, K=k, model="aspectRcpp", method="VB",
                                  alpha=alpha, beta=beta, gamma=gamma,
                                  iter=vb.iters)
    E_Z <- modelVBaspect$E_Z
    E_W <- modelVBaspect$E_W
    E_H <- modelVBaspect$E_H
    pred <- loglikelihood(modelVBaspect, V.test)
    df.results <- list(date = date(),
                       xphash = hash,
                       xp=xp,
                       dataset = dataset,
                       ntest = ntest,
                       model= "Aspect",
                       K=k,
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
                         xp=xp,
                         dataset = dataset,
                         ntest = ntest,
                         model= "cAspect",
                         K=k,
                         loglikelihood = pred$loglikelihood)
      save_result(file = results_file, df.results)
    }
    
    # NMF-KL -------------------------------------------------------------
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
      df.results <- list(date = date(),
                         xphash = hash,
                         xp=xp,
                         dataset = dataset,
                         ntest = ntest,
                         model= "KL",
                         K=k,
                         loglikelihood = like)
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
                         xp=xp,
                         dataset = dataset,
                         ntest = ntest,
                         model= "logPCA",
                         K=k,
                         loglikelihood = like)
      save_result(file = results_file, df.results)
    }
    
    # LDA ----------------------
    # Warning: Each row of the input matrix needs to contain at least one non-zero entry
    # because otherwise is a non-oberved word
    if(models.LDA){
      V_train_ <- V.train
      V_train_[is.na(V_train_)] <- 0
      # ignore unused features
      idx_unused <- which(apply(V_train_, 1, sum) == 0)
      V_train_ <- V_train_[-idx_unused,]
      res <- LDA(as.simple_triplet_matrix(V_train_), k=k, method = "VEM")
      like <- logLik(res,as.simple_triplet_matrix(V.test))
      df.results <- list(date = date(),
                         xphash = hash,
                         xp=xp,
                         dataset = dataset,
                         ntest = ntest,
                         model= "LDA",
                         K=k,
                         loglikelihood = like)
      save_result(file = results_file, df.results)
    }
    
  }

  # Logistic PCA -----------------------------------------------
  # We do not use this implementation
  if(models.logPCA){
    res <- cv.lsvd(V.train, ks=c(2:9), max_iter=2000) # CV to choose K
    K <- as.numeric(rownames(res))[which.max(res)]
    logpca_model = logisticPCA::logisticPCA(V.train, k=15)
    W <- logpca_model$PCs
  }

}

# Plots ------------------------------------------------------------------------
#

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
df$model[df$model == "KL"]             <- "7.KL"
df$model[df$model == "LDA"]            <- "8.LDA"

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
for(i in 1:nrow(df)){
  if(df$model[i] == "7.KL") {
    df$model[i] <- paste0("7.KL-", df$K[i])
  }
}

for(i in 1:nrow(df)){
  if(df$model[i] == "8.LDA") {
    df$model[i] <- paste0("8.LDA-", df$K[i])
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