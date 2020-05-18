###############################################################################
# Analyse of the influence of K in the non-parametric models
###############################################################################
library(logisticPCA)
library(digest)
library(foreach)
library(doParallel)
library(bigmemory)
library(Matrix)

devtools::load_all()

data("animals")
V <- animals
dataset <- "animals"

gibbs.samples <- 5000
burnin <- 0.9
vb.iters      <- 500
repetitions   <- 3
df.results <- data.frame()
for(k in seq(110,300, by=25)){
  # Gibbs DirBer --------------------------------------------------------------------
  alpha = 1
  beta = 1
  gamma = 1
  modelGibbsDirBer <- BernoulliNMF(V, K=k, model="DirBer", 
                                   alpha=alpha, beta=beta, gamma=gamma/k, 
                                   iter=gibbs.samples, burnin = burnin)
  
  E_Z <- modelGibbsDirBer$E_Z
  E_W <- modelGibbsDirBer$E_W
  E_H <- modelGibbsDirBer$E_H
  pred <- loglikelihood(modelGibbsDirBer, V)
  df.results <- bind_rows(df.results, list(model= "Gibbs Beta-Dir",
                                           loglikelihood = pred$loglikelihood,
                                           K=k,
                                           dataset = dataset))
}

kseq <- seq(1,200, by=10)
kseq <- seq(2,9)
kseq <- seq(200, 1000, by=100)
for(k in kseq){
  # Gibbs DirBer --------------------------------------------------------------------
  alpha = 1
  beta = 1
  gamma = 1
  modelVBDirBer <- BernoulliNMF(V, K=k, model="DirBer", method="VB",
                                   alpha=alpha, beta=beta, gamma=gamma/k, 
                                   iter=vb.iters)
  
  E_Z <- modelVBDirBer$E_Z
  E_W <- modelVBDirBer$E_W
  E_H <- modelVBDirBer$E_H
  pred <- loglikelihood(modelVBDirBer, V)
  df.results <- bind_rows(df.results, list(model= "VB Beta-Dir",
                                           loglikelihood = pred$loglikelihood,
                                           K=k,
                                           dataset = dataset))
}

ggplot(df.results, aes(x=K, y=loglikelihood, color=model)) + geom_point()