# Compute and show the reconstructed matrix 
# (the best, if the algorithm has to choose K)

devtools::load_all()
library(dplyr)
library(tidyr)
library(ggplot2)
library(Matrix)
library(gtools)
#library(data.table)
library(logisticPCA)
library(unvotes)

text_components <-  7
text_features <- 7

data("unvotes100coldwar_absna")
data("paleo")
data("catalanparliament")

# Manually add some absent country codes
un_votes[un_votes$country == 'Federal Republic of Germany', ]$country_code <- 'FG'
un_votes[un_votes$country == 'Zanzibar', ]$country_code <- 'ZN'
un_votes[un_votes$country == 'Yemen Arab Republic', ]$country_code <- 'YR'
#country_codes <- un_votes$country_code[match(rownames(V), un_votes$country)]

# Datasets used in this experiment in the paper
dataset_names <- c('paleo', 'unvotes100coldwar_absna', 'parlamentcat')
datasets <- list(paleo, unvotes100coldwar_absna, catalanparliament)

# Datasets to do a fast check
dataset_names <- c("parlamentcat")
datasets <- list(catalanparliament)

# Training parameters
# Use small parameters for fast checks 
# Use larger parameters for accurate results 
gibbs.samples <- 100
burnin <- 0.9
vb.iters <- 100

# Hyper-parameters
K = 100
alpha = 1
beta = 1
gamma = 1

df.results <- data.frame()
for (i in 1:length(dataset_names)) {
  dataset <- dataset_names[i]
  V <- datasets[[i]]
  F <- dim(V)[1]
  N <- dim(V)[2]
  
  plot_V(V, rowlabels=FALSE)
  
  # Gibbs DirBer ---------------------------------------------------------------
  res <- BernoulliNMF(V, 
                      K=K, 
                      model="DirBer", 
                      alpha=1, 
                      beta=1, 
                      gamma=1/K, 
                      iter=gibbs.samples,
                      burnin = burnin)
  # Expectations
  E_Z <- res$E_Z
  E_W <- res$E_W
  E_H <- res$E_H
  E_V <- expectation(res)
  
  # Estimated log-likelihood from samples of the posterior of V. 
  like <- loglikelihood(res, V)$loglikelihood
  
  # Save results
  filename <- paste0("fig_E_V_GibbsDirBer_reconstruc_", dataset, ".eps")
  p <- plot_V(E_V)
  ggsave(p, 
         filename = filename, 
         height=20, 
         width=20, 
         units='cm')
  df.results <- bind_rows(df.results, list(dataset = dataset, 
                                           model = "Gibbs Beta-Dir", 
                                           likelihood = like))
  
  # Optional: 
  # * A single sample V, W, H free from label switching
  # (Monte-Carlo estimations of E[W], E[H] may contain label switchs). 
  # * Trace of the component sizes in the Monte-Carlo samples
  # * Expectation of the dictionary.
  sample <- sample_VWH.GibbsDirBer(res) 
  plot_trace_size_components(res)
  plot_dictionary(E_W, Kmax=10, rowlabels=TRUE, aspect.ratio = 7)
  
  
  # Gibbs DirDir -----------------------------------------------------------------
  res <- BernoulliNMF(V, 
                      K=K, 
                      model="DirDir", 
                      alpha=1, 
                      gamma=1/K, 
                      iter=gibbs.samples, 
                      burnin = burnin)
  
  # Expectations
  E_Z <- res$E_Z
  E_C <- res$E_C
  E_W <- res$E_W
  E_H <- res$E_H
  E_V <- expectation(res)
  
  # Estimated log-likelihood from samples of the posterior of V. 
  like <- loglikelihood(res, V)$loglikelihood

  # Save results
  df.results <- bind_rows(df.results, list(dataset = dataset, 
                                           model = "Gibbs Dir-Dir", 
                                           likelihood = like))
  filename <- paste0("fig_E_V_GibbsDirDir_reconst_", dataset, ".eps")
  p <- plot_V(E_V)
  ggsave(p, 
         filename = filename, 
         height = 20, 
         width = 20, 
         units = 'cm')

  # Optional: 
  # * A single sample V, W, H free from label switching
  # (Monte-Carlo estimations of E[W], E[H] may contain label switchs). 
  # * Trace of the component sizes in the Monte-Carlo samples
  # * Expectation of the dictionary.
  sample <- sample_VWH.GibbsDirDir(res)
  plot_trace_size_components(res)
  plot_dictionary(E_W, Kmax=15, rowlabels=FALSE)
  
  
  # VB DirBer --------------------------------------------------------------------
  res <- BernoulliNMF(V, 
                      K=K, 
                      model="DirBer", 
                      method="VB",
                      alpha=1, 
                      beta=1, 
                      gamma=1/K, 
                      iter=vb.iters)
  
  # Expectations
  E_Z <- res$E_Z
  E_W <- res$E_W
  E_H <- res$E_H
  E_V <- expectation(res)
  
  # Estimated log-likelihood from the variational posterior of V. 
  like <- loglikelihood(res, V)$loglikelihood
  
  # Save results
  p <- plot_V(E_V)
  filename <-  paste0("fig_dmkd_E_V_VBDirBer_reconstruction_", dataset, ".eps")
  ggsave(p, 
         filename = filename, 
         height=20,
         width=20, 
         units='cm')
  lb   <- max(res$lb)
  df.results <- bind_rows(df.results, list(dataset = dataset, 
                                           model = "VB Beta-Dir", 
                                           likelihood = like,
                                           lowerbound = lb))
  
  # Optional: 
  # * Expectation of the dictionary.
  plot_dictionary(E_W, Kmax=10, rowlabels=FALSE)
  
  
  # VB Aspect --------------------------------------------------------------------
  Kmax <- 9
  likes <- rep(-Inf, Kmax)
  lowerbounds <- rep(-Inf, Kmax)
  models <- list()
  for(k in 1:Kmax){
    res <- BernoulliNMF(V, K=k, model="aspectRcpp", method="VB",
                                  alpha=1, beta=1, gamma=1, 
                                  iter=vb.iters)
    # Expectations
    E_Z <- res$E_Z
    E_W <- res$E_W
    E_H <- res$E_H
    
    # Estimated log-likelihood from samples of the posterior of V. 
    like <- loglikelihood(res, V)$loglikelihood
    
    # Save results for this number of components k.
    lb   <- max(res$lowerbounds)
    likes[k] <- like
    lowerbounds[k] <- lb
    models[[k]] <- res
    df.results <- bind_rows(df.results, list(dataset = dataset, 
                                             model = "VB Aspect", 
                                             likelihood = like,
                                             lowerbound = lb,
                                             K=k))
    
  }
  
  # Choose the best model
  res <- models[[which.max(likes)]]
  
  # Expectations
  E_Z <- res$E_Z
  E_W <- res$E_W
  E_H <- res$E_H
  E_V <- expectation(res)

  # Save results
  filename <- paste0("fig_E_V_VBAspect_reconstruction_", dataset, ".eps")
  p <- plot_V(E_V)
  ggsave(p, 
         filename = filename, 
         height=20, width=20, units='cm')
  
  # Optional: 
  # * Expectation of the dictionary.
  plot_dictionary(E_W, Kmax=100, rowlabels=TRUE)
  
  
  # logPCA --------------------------------------------------------------------
  Kmax <- 20
  likes <- rep(-Inf, Kmax)
  models <- list()
  for(k in 1:Kmax){
    logpca_model <- logisticPCA::logisticSVD(V, k, max_iters = 2000)
    W <- logpca_model$A
    H <- logpca_model$B
    E_V <- fitted(logpca_model, type = "response")
    V.probs <- V*E_V + (1-V)*(1-E_V)
    like <- sum(log(V.probs), na.rm = TRUE)
    likes[k] <- like
    models[[k]] <- logpca_model
    df.results <- bind_rows(df.results, list(dataset = dataset,
                                             model= "logPCA",
                                             loglikelihood = like,
                                             K=k))
  }
  
  # Choose the bets model
  res <- models[[which.max(likes)]]
  
  # Expectation of the posterior of V
  E_V <- fitted(res, type = "response")
  
  # Save results
  p <- plot_V(E_V)
  filename <- paste0("fig_dmkd_E_V_logPCA_reconstruction_", dataset, ".eps")
  ggsave(p, 
         filename = filename, 
         height=20, 
         width=20, 
         units='cm')
}