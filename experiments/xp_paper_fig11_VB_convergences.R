# Note: c-bICA and bICA use the same Z init.
# Beta-Dir VB cannot use it because we initialize it with 
# a large K that approximates infinity (K=100)

devtools::load_all()
library(digest)


data("unvotes100coldwar_absna")
data("paleo")
data("lastfm")
data("catalanparliament")
data("animals")

V <- unvotes100coldwar_absna
dataset <- "unvotes100coldwar_absna"
k <- 4
K <- 100

V <- paleo
dataset <- "paleo"
K <- 100
k <- 4

V <- lastfm
dataset <- "lastfm"
k <- 5
K <- 100

V <- catalanparliament
dataset <- "parliament"
k <- 5
K <- 100

V <- animals
dataset <- "animals"
k <- 5
K <- 100

F <- dim(V)[1]
N <- dim(V)[2]

df.results <- data.frame()

repetitions <- 30
for(xp in 1:repetitions){
  
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
  
  
  for(i in c(2, seq(10,100, by=10))){
  
    # Random initialization of Z matrix
    K_init <- k
    F <- dim(V)[1]
    N <- dim(V)[2]
    Z <- array(NA, dim=c(F, N))
    for(n in 1:N){
      for(f in 1:F){
        if(!is.na(V.train[f,n])){
          Z[f,n] <- sample(K_init,1)
        }
      }
    }
    
    
      
    # Aspect -------------------------------------------------------------------
    modelVBaspect <- BernoulliNMF(V.train, K=k, model="aspectRcpp", method="VB",
                                  alpha=1, beta=1, gamma=1, 
                                  iter=i, Z_init = Z)
    
    pred <- loglikelihood(modelVBaspect, V.test)
    
    df.results <- bind_rows(df.results, list(xphash = hash,
                                             xp=xp, 
                                             model= "Aspect",
                                             loglikelihood = pred$loglikelihood,
                                             K=k,
                                             iter = i,
                                             dataset = dataset,
                                             ntest = ntest))
    
    # collapsed Aspect ---------------------------------------------------------
    modelVBDirBer <- BernoulliNMF(V.train, K=k, model="DirBer", method="VB",
                                  alpha=1, beta=1, gamma=1, 
                                  iter=i, Z_init = Z)
    
    pred <- loglikelihood(modelVBDirBer, V.test)
    df.results <- bind_rows(df.results, list(xphash = hash,
                                             xp=xp, 
                                             model= "collapsed Aspect",
                                             loglikelihood = pred$loglikelihood,
                                             K=k,
                                             iter = i,
                                             dataset = dataset,
                                             ntest = ntest))
    
    
    # Beta-Dir init k small ----------------------------------------------------
    modelVBDirBer <- BernoulliNMF(V.train, K=K, model="DirBer", method="VB",
                                  K_init = k,
                                  alpha=1, beta=1, gamma=1/K, 
                                  iter=i, Z_init = Z)
    
    pred <- loglikelihood(modelVBDirBer, V.test)
    df.results <- bind_rows(df.results, list(xphash = hash,
                                             xp=xp, 
                                             model= "VB Beta-Dir - k small",
                                             loglikelihood = pred$loglikelihood,
                                             K=K,
                                             iter = i,
                                             dataset = dataset,
                                             ntest = ntest))
   
    # Beta-Dir ---------------------------------------------------------
    modelVBDirBer <- BernoulliNMF(V.train, K=K, model="DirBer", method="VB",
                                  alpha=1, beta=1, gamma=1/K, 
                                  iter=i)
    
    pred <- loglikelihood(modelVBDirBer, V.test)
    df.results <- bind_rows(df.results, list(xphash = hash,
                                             xp=xp, 
                                             model= "VB Beta-Dir",
                                             loglikelihood = pred$loglikelihood,
                                             K=K,
                                             iter = i,
                                             dataset = dataset,
                                             ntest = ntest))
    
     
  }
}

# Write results to file
write.table(df.results, file = "results_dmkd_vb_convergence.csv",
            append = TRUE, col.names = FALSE, 
            quote = FALSE, sep = ",", row.names = FALSE)



# Plot ------------------------------------------------------------------------
df <- read.table( file = "results_dmkd_vb_convergence.csv", sep=",")
names(df) <- c('xphash', 'xp', 'model', 'loglikelihood', 
               'k', 'iter', 'dataset', 'ntest')

df <- df %>% mutate(perplexity = -loglikelihood/ntest)
df <- df %>% group_by(dataset) %>% 
  mutate(rel_loglikelihood = loglikelihood - max(loglikelihood))

df$model <- as.character(df$model)
df$dataset <- as.character(df$dataset)
df$dataset[grepl("unvotes", df$dataset)] <- "unvotes"
df$model[df$model == "Aspect"] <- "bICA"
df$model[grepl("collapsed Aspect", df$model)] <- "c-bICA"
df$model[df$model == "VB Beta-Dir"] <- "Beta-Dir VB"
df$model[grepl("k small", df$model)] <- "Beta-Dir VB*"
df$model   <- as.factor(df$model)
df$dataset <- as.factor(df$dataset)

df <- df %>% filter(model != "Beta-Dir VB*")
# Boxplots
base_size <- 8
p <- ggplot(df, aes(x=as.factor(iter), y=perplexity)) + 
  geom_boxplot(outlier.shape = NA, aes(colour=model)) + 
  facet_grid(dataset ~ ., scales = "free") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, size=base_size*1, hjust = 1, colour = "black"),
        axis.text.y = element_text(size= base_size, color = 'black'),
        strip.background =element_rect(fill="white"),
        aspect.ratio = 1/3) + xlab('iterations')
#plot.margin=unit(c(0,1,0,-3), "mm")
print(p)

# Mean lines
base_size <- 8
df <- df %>% group_by(dataset, model, iter) %>% 
  summarise(mean = mean(perplexity), sd = sd(perplexity))
p <- ggplot(df %>% filter(iter <100), 
            aes(x=iter, y=mean, color = model, shape= model, linetype=model)) + 
  geom_point() +
  geom_line() +
  facet_grid(dataset ~ ., scales = "free") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 0, size=base_size*1, hjust = 1, colour = "black"),
        axis.text.y = element_text(size= base_size, color = 'black'),
        strip.background =element_rect(fill="white"),
        legend.title = element_blank(),
        legend.position = "top",
        #legend.text=element_text(size=3),
        #legend.position = c(0.85, 0.1),
        #legend.background = element_rect(color = "black"),
        aspect.ratio = 1/3.5) + xlab('iterations') + ylab('perplexity')
#plot.margin=unit(c(0,1,0,-3), "mm")
print(p)

ggsave(p, filename = "fig_dmkd_vb_convergence.eps", height=17, width=13, units='cm')


