# For many datasets and many algorithms
# show the dictionary of each algorithm in each dataset
# (the best, if the algorithm has to choose K)

devtools::load_all()
library(dplyr)
library(tidyr)
library(ggplot2)
library(Matrix)
library(gtools)
library(logisticPCA)
library(unvotes)

text_components <-  7
text_features <- 7

data("animals")
data("paleo")
data("catalanparliament")
data(catalanparliament_labels) # we use this to add the party to each MP.
data("unvotes100coldwar_absna")

# Datasets used in this experiment in the paper
dataset_names <- c('paleo', 'unvotes100coldwar_absna', 'parlamentcat')
datasets <- list(paleo, unvotes100coldwar_absna, catalanparliament)

# Datasets to do a fast check
dataset_names <- c("animals", "parlamentcat")
datasets <- list(animals, catalanparliament)

# Training parameters
# Use small parameters for fast checks 
# Use larger parameters for accurate results 
gibbs.samples <- 10
burnin <- 0.9
vb.iters <- 10

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
  
  # Gibbs DirBer -----------------------------------------------------------------
  modelGibbsDirBer <- BernoulliNMF(V, 
                                   K=K, 
                                   model="DirBer", 
                                   alpha=alpha, 
                                   beta=beta, 
                                   gamma=gamma/K, 
                                   iter=gibbs.samples, 
                                   burnin = burnin)
  E_W <- modelGibbsDirBer$E_W
  like <- loglikelihood(modelGibbsDirBer, V)$loglikelihood # quality of fit
  df.results <- bind_rows(df.results, list(dataset = dataset, 
                                           model = "Gibbs Beta-Dir", 
                                           likelihood = like))
  p <- plot_dictionary(E_W, Kmax=10, rowlabels=TRUE, aspect.ratio = 7)

  # Gibbs DirDir -----------------------------------------------------------------
  modelGibbsDirDir <- BernoulliNMF(V, 
                                   K=K, 
                                   model="DirDir", 
                                   alpha=alpha,
                                   beta=beta,
                                   gamma=gamma/K, 
                                   iter=gibbs.samples, 
                                   burnin = burnin)
  
  E_W <- modelGibbsDirDir$E_W
  like <- loglikelihood(modelGibbsDirDir, V)$loglikelihood
  df.results <- bind_rows(df.results, list(dataset = dataset, 
                                           model = "Gibbs Dir-Dir", 
                                           likelihood = like))
  p <- plot_dictionary(E_W, Kmax=15, rowlabels=FALSE)

  # VB Aspect --------------------------------------------------------------------
  Kmax <- 15
  likes <- rep(-Inf, Kmax)
  lowerbounds <- rep(-Inf, Kmax)
  models <- list()
  for(k in 1:Kmax){
    modelVBaspect <- BernoulliNMF(V, 
                                  K=k, 
                                  model="aspectRcpp", 
                                  method="VB",
                                  alpha=alpha, 
                                  beta=beta, 
                                  gamma=gamma, 
                                  iter=vb.iters)
    E_W <- modelVBaspect$E_W
    like <- loglikelihood(modelVBaspect, V)$loglikelihood
    lb   <- max(modelVBaspect$lowerbounds)
    likes[k] <- like
    lowerbounds[k] <- lb
    models[[k]] <- modelVBaspect
    df.results <- bind_rows(df.results, list(dataset = dataset, 
                                             model = "VB Aspect", 
                                             likelihood = like,
                                             lowerbound = lb,
                                             K=k))
    
  }
  modelVBaspect <- models[[which.max(likes)]]
  E_W <- modelVBaspect$E_W
  p <- plot_dictionary(E_W, Kmax=100, rowlabels=TRUE)

  
  # logPCA --------------------------------------------------------------------
  Kmax <- 8 # set it to K in the aspect model
  likes <- rep(-Inf, Kmax)
  models <- list()
  for(k in 1:Kmax){
    logpca_model <- logisticPCA::logisticSVD(V, k, max_iters = 2000)
    W <- logpca_model$A
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
  modelLogPCA <- models[[which.max(likes)]]
  W <- modelLogPCA$A
  p <- plot_dictionary(W, Kmax=100, rowlabels=TRUE)
  
  
  # VB DirBer --------------------------------------------------------------------
  modelVBDirBer <- BernoulliNMF(V, 
                                K=K, 
                                model="DirBer", 
                                method="VB",
                                alpha=1, 
                                beta=1, 
                                gamma=1/K, 
                                iter=vb.iters)
  
  E_W <- modelVBDirBer$E_W
  like <- loglikelihood(modelVBDirBer, V)$loglikelihood
  lb   <- max(modelVBDirBer$lb)
  df.results <- bind_rows(df.results, list(dataset = dataset, 
                                           model = "VB Beta-Dir", 
                                           likelihood = like,
                                           lowerbound = lb))
  p <- plot_dictionary(E_W, Kmax=10, rowlabels=FALSE)

  
  # Plots ----------------------------------------------------------------------
  # Each dictionary neeeds different aesthetics
   if (dataset == 'animals'){
     labels <- rownames(modelGibbsDirBer$V)
     idx.active <- sort(order(-rowSums(modelGibbsDirBer$V))[1:100]) # by freq.
     parameters = list(idx.subset = idx.active,
                       labels = labels,
                       Kmax = 8,
                       aspect.ratio = 10.5,
                       ggsave.height=34, 
                       ggsave.width=20)
   } else if (dataset == 'parlamentcat') {
     parties <- catalanparliament_labels$party[match(rownames(modelGibbsDirBer$V), 
                                                     catalanparliament_labels$MP)]
     levels(parties) <- c("CEC", "Cs", "CUP", "ERC", "JxC", "PP", "PSC")
     labels <- paste(rownames(modelGibbsDirBer$V), '-', parties)
     idx.active <- sort(order(-rowSums(modelGibbsDirBer$V))[1:300])
     parameters = list(idx.subset = idx.active,
                       labels = labels,
                       Kmax = 8,
                       aspect.ratio = 11,
                       ggsave.height=34, 
                       ggsave.width=20) 
   } else if (dataset == 'paleo') {
     labels <- rownames(modelGibbsDirBer$V)
     idx.active <- sort(order(-rowSums(modelGibbsDirBer$V))[1:100])
     parameters = list(idx.subset = idx.active,
                       labels = SAME,
                       Kmax = 8,
                       aspect.ratio = 10.5,
                       ggsave.height=34, 
                       ggsave.width=20) 
   } else if (dataset == 'unvotes100coldwar_absna') {
     country_shortnames <- rownames(V)
     country_shortnames[country_shortnames == 'United States of America'] <- "U.S.A"
     country_shortnames[grepl("Bolivia", country_shortnames)] <- "Bolivia"
     country_shortnames[grepl("United Kingdom", country_shortnames)] <- "United Kingdom"
     country_shortnames[grepl("Venezuela", country_shortnames)] <- "Venezuela"
     country_shortnames[grepl("Iran", country_shortnames)] <- "Iran"
     country_shortnames[grepl("Syrian", country_shortnames)] <- "Syria"
     country_shortnames[grepl("Tanzania", country_shortnames)] <- "Tanzania"
     country_shortnames[grepl("Lao", country_shortnames)] <- "Laos"
     country_shortnames[grepl("Congo", country_shortnames)] <- "Congo"
     country_shortnames[grepl("Yemen Arab Republic", country_shortnames)] <- "Yemen Aran Rep."
     country_shortnames[grepl("Russian Federation", country_shortnames)] <- "Russian Fed."
     labels = country_shortnames
     parameters = list(idx.subset = idx.active,
                       labels = labels,
                       Kmax = 8,
                       aspect.ratio = 10.5,
                       ggsave.height=35, 
                       ggsave.width=20) 
   }
  
   # cast logPCA dictionary from -inf,+inf to [0,1] for visualization purposes
   rownames(modelLogPCA$A) <- rownames(modelGibbsDirBer$V)
   A <- (modelLogPCA$A-min(modelLogPCA$A)) / (max(modelLogPCA$A)-min(modelLogPCA$A))
   
   p <- plot_multiple_dictionaries(list(modelGibbsDirBer$E_W,
                                        modelVBDirBer$E_W,
                                        modelGibbsDirDir$E_W,
                                        modelVBaspect$E_W,
                                        A),
                                   names = c("Beta-Dir GS",
                                             "Beta-Dir VB",
                                             "Dir-Dir GS",
                                             "bICA",
                                             "logPCA"),
                                   labels = parameters$labels,
                                   rowlabels = TRUE,
                                   sort = TRUE,
                                   Kmax = parameters$Kmax,
                                   aspect.ratio = parameters$aspect.ratio,
                                   idx.subset = parameters$idx.active)
  ggsave(p, 
         filename = paste0("fig_dictionaries_", dataset, ".eps"), 
         height=parameters$ggsave.height, 
         width=parameters$ggsave.width, 
         units='cm')
}