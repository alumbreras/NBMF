# The authors should investigate visualizing the latent space in 2D 
# by restricting to two components or by applying PCA. 
# The resulting scatter plots can provide further interpretable insides 
# and demonstrate advantages of binary PCA.

devtools::load_all()
library(dplyr)
library(tidyr)
library(ggplot2)
library(Matrix)
library(gtools)
library(data.table)
library(logisticPCA)

library(unvotes)
data("unvotes100coldwar_absna")
V <- unvotes100coldwar_absna
un_votes[un_votes$country == 'Federal Republic of Germany', ]$country_code <- 'FG'
un_votes[un_votes$country == 'Zanzibar', ]$country_code <- 'ZN'
un_votes[un_votes$country == 'Yemen Arab Republic', ]$country_code <- 'YR'
country_codes <- un_votes$country_code[match(rownames(V), un_votes$country)]
dataset <- 'unvotes100coldwar_absna'

V <- t(V) # I want H to represent countries-topics with Beta ranges [0,1]
names <- country_codes

modelVBDirBer <- BernoulliNMF(V, K=1, model="DirBer", method="VB",
                              alpha=1, beta=1, gamma=1/2, 
                              iter=1000)
E_W <- modelVBDirBer$E_W
E_H <- modelVBDirBer$E_H

logpca_model <- logisticPCA::logisticSVD(V, k=2, max_iters = 2000)
W <- logpca_model$A
H <- logpca_model$B

par(mfrow=c(1,1))
plot(E_H, cex= 0, main="H coefficients in Beta-Dir")
text(E_H[,1], E_H[,2], labels=names, cex= 0.7)

plot(H, cex= 0, main="H coefficients in log-PCA")
text(H[,1], H[,2], labels=names, cex= 0.7)

# Plot a logPCA dictionary with countries-topics to show that plotting in the log space
# is useless
logpca_model <- logisticPCA::logisticSVD(V, k=8, max_iters = 2000)
W <- logpca_model$A
H <- logpca_model$B

A_mapped <- (logpca_model$A-min(logpca_model$A)) / (max(logpca_model$A)-min(logpca_model$A))

p <- plot_dictionary(logpca_model$A, Kmax=8, rowlabels=TRUE)
p <- plot_dictionary(A_mapped, Kmax=8, rowlabels=TRUE)


