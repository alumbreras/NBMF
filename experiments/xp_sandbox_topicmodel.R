library(topicmodels) # LDA (TODO install it and use perplexity function)
library(tm)
library(slam)
results_file = "results_predictions_criteo.csv"
devtools::load_all()

data("catalanparliament")
V <- catalanparliament
dataset <- "parlament"

# See here to create a Document Term Matrix
df.V <- as.data.frame(V)
df.V$from <- rownames(V)
df.V <- gather(df.V, to, count, -from)
df.V.train <- 
df.V.test  <-  
# with triplets
V <- as.simple_triplet_matrix(V)


#https://cran.r-project.org/web/packages/tidytext/vignettes/tidying_casting.html
#https://stackoverflow.com/questions/21355156/topic-models-cross-validation-with-loglikelihood-or-perplexity
res <- LDA(V, k=10, method = "VEM")
perplexity(res, V)
logLik(res,V)

full_data  <- AssociatedPress
n <- nrow(full_data)
#-----------validation--------
k <- 5

splitter <- sample(1:n, round(n * 0.75))
train_set <- full_data[splitter, ]
valid_set <- full_data[-splitter, ]


lda_inf <- posterior(lda, AssociatedPress[21:30,])
# perplexity