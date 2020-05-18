# Components size --------------------------------------------------------------
plot_trace_Z_norms <- function(Z_samples){
  # Plot trace of Z norms (size of every dimension)
  J <- dim(Z_samples)[3]
  df.components_size <- NULL #data.frame(iteration, component, size)
  Kmax <- max(Z_samples)
  for(j in 1:J){
    cat('\niter:', j)
    freqs <- table(Z_samples[,,j])
    for(k in 1:length(freqs)){
      size <- freqs[k]
      row <- data.frame(iteration=j, component=k, size=size)
      df.components_size <- rbindlist(list(df.components_size, row))
    }
  }
  
  
  p <- ggplot(df.components_size, aes(x=iteration, y=size, group=component)) + 
    geom_line() +
    theme_bw() +
    theme(legend.position = "none",
          aspect.ratio=1) +   coord_fixed() 
  print(p)
  ggsave(p, filename = "Z_sizes_parlamentcat.eps", height=8, width=8, units='cm')
}


#TODO
# Plot trace of number of components used
# J <- dim(samples$Z_samples)[3]
# df.n_components <- c(NA, NA)
# for(j in 1:J){
#   df.n_components <- rbind(df.n_components, c(j, length(table(samples$Z_samples[,,1])))) 
# }
