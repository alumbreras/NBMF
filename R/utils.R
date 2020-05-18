save_result <- function(filename, results){
  write.table(results, file = filename,
              append = TRUE, col.names = !file.exists(filename),
              quote = FALSE, sep = ",", row.names = FALSE)
}