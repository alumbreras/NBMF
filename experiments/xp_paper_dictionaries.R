# For many datasets and many algorithms
# show the dictionary of each algorithm in each dataset
# (the best, if the algorithm has to choose K)
################################################################################
devtools::load_all()
library(dplyr)
library(tidyr)
library(ggplot2)
library(Matrix)
library(gtools)
library(data.table)
library(logisticPCA)

text_components <-  7
text_features <- 7

data("animals")
V <- animals
dataset <- 'animals'



data("lastfm")
V <- lastfm
dataset <- "lastfm"

library(unvotes)
data("unvotes100coldwar_absna")
V <- unvotes100coldwar_absna
un_votes[un_votes$country == 'Federal Republic of Germany', ]$country_code <- 'FG'
un_votes[un_votes$country == 'Zanzibar', ]$country_code <- 'ZN'
un_votes[un_votes$country == 'Yemen Arab Republic', ]$country_code <- 'YR'
country_codes <- un_votes$country_code[match(rownames(V), un_votes$country)]
dataset <- 'unvotes100coldwar_absna'


data("paleo")
V <- paleo
dataset <- 'paleo'

data("catalanparliament")
data(catalanparliament_labels)
V <- catalanparliament
dataset <- 'parlamentcat'

F <- dim(V)[1]
N <- dim(V)[2]

# Plot V
plot_V(V, rowlabels=TRUE, collabels = FALSE)
p <- plot_V(V, rowlabels=FALSE)

K = 100
alpha = 1
beta = 1
gamma = 1

gibbs.samples <- 3000
burnin <- 0.9
vb.iters <- 1000

df.results <- data.frame()

# Gibbs DirBer -----------------------------------------------------------------
modelGibbsDirBer <- BernoulliNMF(V, K=K, model="DirBer", 
                                 alpha=1, beta=1, gamma=1/K, 
                                 iter=gibbs.samples, burnin = burnin)
E_Z <- modelGibbsDirBer$E_Z
E_W <- modelGibbsDirBer$E_W
E_H <- modelGibbsDirBer$E_H
sample <- sample_VWH.GibbsDirBer(modelGibbsDirBer) # sample W, H to avoid label switching of E[W]
like <- loglikelihood(modelGibbsDirBer, V)$loglikelihood # quality of fit
df.results <- bind_rows(df.results, list(dataset = dataset, 
                                         model = "Gibbs Beta-Dir", 
                                         likelihood = like))

# Matrix E[W]
p <- plot_dictionary(E_W, Kmax=10, rowlabels=TRUE, aspect.ratio = 7)
# Reconstructed E[V]
E_V <- expectation(modelGibbsDirBer)
p <- plot_V(E_V)
#plot_trace_size_components(modelGibbsDirBer)

# Gibbs DirDir -----------------------------------------------------------------
modelGibbsDirDir <- BernoulliNMF(V, K=K, model="DirDir", 
                                 alpha=1, gamma=1/K, 
                                 iter=gibbs.samples, burnin = burnin)



E_Z <- modelGibbsDirDir$E_Z
E_C <- modelGibbsDirDir$E_C
E_W <- modelGibbsDirDir$E_W
E_H <- modelGibbsDirDir$E_H
sample <- sample_VWH.GibbsDirBer(modelGibbsDirDir)
plot_trace_size_components(modelGibbsDirDir)

like <- loglikelihood(modelGibbsDirDir, V)$loglikelihood # quality of fit
df.results <- bind_rows(df.results, list(dataset = dataset, 
                                         model = "Gibbs Dir-Dir", 
                                         likelihood = like))

# Matrix E[W]
p <- plot_dictionary(E_W, Kmax=15, rowlabels=FALSE)
# Reconstructed E[V]
E_V <- expectation(modelGibbsDirDir)
p <- plot_V(E_V)

# VB Aspect --------------------------------------------------------------------
Kmax <- 15
likes <- rep(-Inf, Kmax)
lowerbounds <- rep(-Inf, Kmax)
models <- list()
for(k in 1:Kmax){
  modelVBaspect <- BernoulliNMF(V, K=k, model="aspectRcpp", method="VB",
                                alpha=1, beta=1, gamma=1, 
                                iter=vb.iters)
  E_Z <- modelVBaspect$E_Z
  E_W <- modelVBaspect$E_W
  E_H <- modelVBaspect$E_H
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
E_Z <- modelVBaspect$E_Z
E_W <- modelVBaspect$E_W
E_H <- modelVBaspect$E_H

# Squared W with E[W]
p <- plot_dictionary(E_W, Kmax=100, rowlabels=TRUE)
E_V <- expectation(modelVBaspect)
p <- plot_V(E_V)


# logPCA --------------------------------------------------------------------
Kmax <- 8 # set it to K in the aspect model
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

modelLogPCA <- models[[which.max(likes)]]
W <- modelLogPCA$A
H <- modelLogPCA$B

# Squared W with E[W]
p <- plot_dictionary(W, Kmax=100, rowlabels=TRUE)


# VB DirBer --------------------------------------------------------------------
modelVBDirBer <- BernoulliNMF(V, K=K, model="DirBer", method="VB",
                              alpha=1, beta=1, gamma=1/K, 
                              iter=vb.iters)

E_Z <- modelVBDirBer$E_Z
E_W <- modelVBDirBer$E_W
E_H <- modelVBDirBer$E_H
like <- loglikelihood(modelVBDirBer, V)$loglikelihood
lb   <- max(modelVBDirBer$lb)
df.results <- bind_rows(df.results, list(dataset = dataset, 
                                         model = "VB Beta-Dir", 
                                         likelihood = like,
                                         lowerbound = lb))


# Squared W with E[W]
p <- plot_dictionary(E_W, Kmax=10, rowlabels=FALSE)
E_V <- expectation(modelVBDirBer)
p <- plot_V(E_V)


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


# Plot ONU dictionaries --------
# dictionary -inf,+inf -> [0,1]
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
                                rowlabels = TRUE,
                                sort = TRUE,
                                Kmax=8,
                                aspect.ratio = 10.5, labels = country_shortnames)
ggsave(p, filename = paste0("fig_dmkd_dictionaries_", dataset, ".eps"), 
       height=35, width=20, units='cm')


# Plot paleo ----------------
idx.active <- sort(order(-rowSums(modelGibbsDirBer$V))[1:100])# by freq.

# dictionary -inf,+inf -> [0,1]
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
                                rowlabels = TRUE,
                                sort = TRUE,
                                Kmax=8,
                                aspect.ratio = 10.5, idx.subset = idx.active)
ggsave(p, filename = paste0("fig_dmkd_dictionaries_", dataset, ".eps"), 
       height=34, width=20, units='cm')


# Plot parliament --------------
idx.active <- sort(order(-rowSums(modelGibbsDirBer$V))[1:100])# by freq.

# dictionary -inf,+inf -> [0,1]
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
                                rowlabels = TRUE,
                                sort = TRUE,
                                Kmax=8,
                                aspect.ratio = 10.5, 
                                idx.subset = idx.active)
ggsave(p, filename = paste0("fig_dmkd_dictionaries_", dataset, ".eps"), 
       height=34, width=20, units='cm')

# MP - party -------
parties <- catalanparliament_labels$party[match(rownames(modelGibbsDirBer$V), catalanparliament_labels$MP)]
levels(parties) <- c("CEC", "Cs", "CUP", "ERC", "JxC", "PP", "PSC")
labels <- paste(rownames(modelGibbsDirBer$V), '-', parties)
idx.active <- sort(order(-rowSums(modelGibbsDirBer$V))[1:300])# by freq.

## dictionary -inf,+inf -> [0,1]
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
                                labels = labels,
                                rowlabels = TRUE,
                                sort = TRUE,
                                Kmax=8,
                                aspect.ratio = 11, 
                                idx.subset = idx.active)

ggsave(p, filename = paste0("Fig7_dmkd_dictionaries_parties_", dataset, ".eps"), 
       height=34, width=20, units='cm')