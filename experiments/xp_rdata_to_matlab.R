library(tidyr)
library(reshape2)
library('R.matlab')
devtools::load_all()

data("catalanparliament")
V <- catalanparliament
dataset <- "parlament"
rownames(V) <- 1:dim(V)[1]
colnames(V) <- 1:dim(V)[2]
df.data <- melt(V, varnames=c("f", "n"))

# Connvert to Matlab matrix
writeMat("parliament.mat", parliament = V)





# OLD
#################################################
write.table(df.data, file = "parliament.data",
            append = FALSE, col.names = FALSE,
            quote = FALSE, sep = ",", row.names = FALSE)

labels <- 1:length(unique(df[,1]))
write.table(labels, file = "parliament.label",
            append = FALSE, col.names = FALSE,
            quote = FALSE, sep = ",", row.names = FALSE)


length(unique(df.data$f)) == length(labels)

# length(unique(train.data$V1)) = length(train.label$V1)
train.data <- read.table(file = "train.data")
train.label <- read.table(file = "train.label")

unique(train.data)


test.data <- read.table(file = "test.data")
test.label <- read.table(file = "test.label")