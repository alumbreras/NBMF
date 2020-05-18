library(dplyr)
library(tidyr)
library(ggplot2)
library(Matrix)
library(gtools)
library(data.table)
library(svs) # Classic NMFs
library(NMF) # Classic NMFs
library(logisticPCA)
library(topicmodels) # LDA

text_components <-  7
text_features <- 7

data("catalanparliament")
V <- catalanparliament

# Logistic PCA
logpca_model = logisticPCA::logisticPCA(V, k=10)
W <- logpca_model$PCs
plot_dictionary(W, labels = TRUE)

# NMF package, KL
NMF_model <- NMF::nmf(V, 10, "KL")
W <- basis(NMF_model)
H <- coef(NMF_model)
V.hat <- fitted(NMF_model)
plot_dictionary(W, rowlabels = TRUE)
plot_V(V.hat)

# svs package, KL
NMF_model <- svs::fast_nmf_KL(V, k=10, tol = 1e-08)
H <- NMF_model$pos1
W <- t(NMF_model$pos2)
plot_dictionary(W, rowlabels = TRUE)


# Other NMFs
NMF_model <- svs::fast_nmf(V, k=10, tol = 1e-08)
NMF_model <- svs::fast_nmf_Al(V, k=10, tol = 1e-08)
NMF_model <- svs::fast_nmf_Fr(V, k=10, tol = 1e-08)
H <- NMF_model$pos1
W <- t(NMF_model$pos2)
plot_dictionary(W, rowlabels = TRUE)

# pLSA
pLSA_model <- fast_plsa(V, k=10, tol = 1e-08)
W <- pLSA_model$prob2

# LDA
LDA_model <- LDA(V, k=10, method = "VEM")
W <- t(LDA_model@beta)
W <- LDA_model@gamma

# Plot W ======================================================================
text_components <-  7
text_features <- 7

W <- W[,order(-colSums(W, na.rm = TRUE))] # sort columns by norm.
rownames(W) <- rownames(V)
colnames(W) <- 1:dim(W)[2]

df.W <- data.frame(W)
colnames(df.W) <- 1:dim(W)[2]

# numbered names
df.W$name <- 1:dim(W)[1]
df.W$name <- factor(df.W$name, levels = df.W$name)

# text names
df.W$name <- rownames(W)
df.W$name <- factor(df.W$name, levels = df.W$name)

df.W <- gather(df.W, key='dimension', value='value' , -name)
df.W$dimension <- as.numeric(df.W$dimension)


Kmax <- 100
p <- ggplot(df.W %>% filter(dimension <Kmax), aes(x=factor(dimension), y=name)) + 
  geom_tile(aes(fill = value), colour = "white") + 
  scale_fill_gradient(low = "white", high = "steelblue", na.value = 'red') +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  xlab("component") +
  ylab("from") +
  theme_bw() +
  theme(#axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    #axis.ticks.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    aspect.ratio=1) +   coord_fixed() 
print(p)




V_pred <- fitted(logpca_model, type = "response")
# Plot predicted V ============================================================
df.V <- data.frame(V_pred)
colnames(df.V) <- 1:dim(V)[2]

# numbered names
df.V$from <- 1:dim(V)[1]
df.V$from <- factor(df.V$from, levels = df.V$from)


df.V <- gather(df.V, key='to', value='value' , -from)
df.V$to <- as.numeric(df.V$to)



p <- ggplot(df.V, aes(x=to, y=from)) + 
  geom_tile(aes(fill = value), colour = "white") + 
  scale_fill_gradient(low = "white", high = "steelblue", na.value = 'red') +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        aspect.ratio=1) +   coord_fixed() 

print(p)
