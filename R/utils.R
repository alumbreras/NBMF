save_result <- function(filename, results){
  write.table(results, file = filename,
              append = TRUE, col.names = !file.exists(filename),
              quote = FALSE, sep = ",", row.names = FALSE)
}

sample_V <- function(W,H){
  F <- dim(W)[1]
  N <- dim(H)[2]
  V <- array(NA, dim = c(F,N))
  for(f in 1:F){
    for(n in 1:N){
      V[f,n] <- rbinom(1, 1, prob = W[f,]%*%H[,n])
    }
  }
  return(V)
}