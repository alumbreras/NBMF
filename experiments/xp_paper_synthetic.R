#
# Create synthetis data for Beta-Dir and Dir-Dir and learn the dictionaries
# Theory: BetaDir is similar or better than DirDir fitting the data, since it may use
# all dependencies. However, DirDir may generalize better if the data is DirDir
devtools::load_all()
library(gtools)

F <- 100
N <- 100
K <- 4
W <- array(NA, dim = c(F,K))
H <- array(NA, dim = c(K,N))
V.DirDir  <- array(NA, dim = c(F,N))
V.BetaDir <- array(NA, dim = c(F,N))
V.DirBeta <- array(NA, dim = c(F,N))

# Dir-Dir
alpha <- 1
beta <- 1
eta <- 1
for(f in 1:F){
  W[f,] <- rdirichlet(1, rep(alpha,K))
}

for(n in 1:N){
  H[,n] <- rdirichlet(1, rep(eta,K))
}


for(f in 1:F){
  for(n in 1:N){
    V.DirDir[f,n] <- rbinom(1, 1, prob = W[f,]%*%H[,n])
  }
}

plot_V(V.DirDir)
res <- heatmap(V.DirDir)
idx.row <- res$rowInd
idx.col <- res$colInd
V.DirDir <- V.DirDir[idx.row, idx.col]
p <-plot_V(V.DirDir)
dataset <- "toy"

ggsave(p, filename = paste0("fig_V_synthetic_DirDir_nostruct.eps"), height=5, width=5, units='cm')

# Beta-Dir
alpha <- 0.1
beta <- 0.1
eta <- 0.1
for(f in 1:F){
  W[f,] <- rdirichlet(1, rep(alpha,K))
}

for(n in 1:N){
  H[,n] <- rbeta(K, alpha,beta)
}


for(f in 1:F){
  for(n in 1:N){
    V.BetaDir[f,n] <- rbinom(1, 1, prob = W[f,]%*%H[,n])
  }
}

plot_V(V.BetaDir)
res <- heatmap(V.BetaDir)
idx.row <- res$rowInd
idx.col <- res$colInd
V.BetaDir <- V.BetaDir[idx.row, idx.col]
p <- plot_V(V.BetaDir)
dataset <- "toy"
ggsave(p, filename = paste0("fig_V_synthetic_BetaDir_nostruct.eps"), height=5, width=5, units='cm')


# Beta-Dir
alpha <- 0.1
beta <- 0.1
eta <- 0.1
for(f in 1:F){
  W[f,] <- rbeta(K, alpha,beta)
}

for(n in 1:N){
  H[,n] <- rdirichlet(1, rep(alpha,K))
}


for(f in 1:F){
  for(n in 1:N){
    V.DirBeta[f,n] <- rbinom(1, 1, prob = W[f,]%*%H[,n])
  }
}

plot_V(V.DirBeta)
res <- heatmap(V.DirBeta)
idx.row <- res$rowInd
idx.col <- res$colInd
V.DirBeta <- V.DirBeta[idx.row, idx.col]
p <- plot_V(V.DirBeta)
dataset <- "toy"
ggsave(p, filename = paste0("fig_V_synthetic_DirBeta_nostruct.eps"), height=5, width=5, units='cm')


###############################################################################

df.results <- data.frame()

K <-  100
modelGibbsDirBer <- BernoulliNMF(V, K=K, model="DirBer", 
                                 alpha=1, beta=1, gamma=1/K, 
                                 iter=1000, burnin = 0.9)
E_Z <- modelGibbsDirBer$E_Z
E_W <- modelGibbsDirBer$E_W
E_H <- modelGibbsDirBer$E_H
plot_dictionary(E_W, Kmax=10, rowlabels=TRUE, aspect.ratio = 7)
E_V_DirBer <- expectation(modelGibbsDirBer)
plot_V(E_V_DirBer)

like <- loglikelihood(modelGibbsDirBer, V)$loglikelihood # quality of fit
df.results <- bind_rows(df.results, list(dataset = dataset, 
                                         model = "Gibbs Beta-Dir", 
                                         likelihood = like))


K <-  100
modelGibbsDirDir <- BernoulliNMF(V, K=K, model="DirBer", 
                                 alpha=1, gamma=1/K, 
                                 iter=2000, burnin = 0.9)
E_Z <- modelGibbsDirDir$E_Z
E_W <- modelGibbsDirDir$E_W
E_H <- modelGibbsDirDir$E_H
plot_dictionary(E_W, Kmax=10, rowlabels=TRUE, aspect.ratio = 7)
E_V_DirDir <- expectation(modelGibbsDirDir)
plot_V(E_V_DirDir)
like <- loglikelihood(modelGibbsDirDir, V)$loglikelihood # quality of fit
df.results <- bind_rows(df.results, list(dataset = dataset, 
                                         model = "Gibbs Dir-Dir", 
                                         likelihood = like))