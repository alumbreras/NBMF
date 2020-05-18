# Functions to plot data matrix V, dictionary W, activations matrix H, etc
# ############################################################################
library(ggplot2)
library(tidyr)
library(dplyr)

plot_V_discrete <- function(V, rowlabels = FALSE, collabels=FALSE){
  
  base_size <- 7
  axis.text.y <- ggplot2::element_text(size=base_size*1, hjust = 1, 
                                       colour = "black")
  axis.text.x <- ggplot2::element_text(angle = 45, size=base_size*1, hjust = 1, 
                                       colour = "black")
  
  if(is.null(rownames(V))){
    rownames(V) = 1:dim(V)[1]
  }
  
  if (rowlabels == FALSE){ 
    axis.text.y  <- ggplot2::element_blank()
    axis.ticks.y <- ggplot2::element_blank()
  }
  
  if (collabels == FALSE){ 
    axis.text.x  <- ggplot2::element_blank() 
    axis.ticks.x <- ggplot2::element_blank()
  }
  
  
  df.V <- data.frame(V)
  
  # Labs as factors preserving original order
  df.V$rowlab <- rownames(V)
  df.V$rowlab <- factor(df.V$rowlab, levels = df.V$rowlab)
  mytheme <- ggplot2::theme(
    axis.text.x   = axis.text.x,
    axis.text.y   = axis.text.y,
    axis.ticks.x  = element_blank(),
    axis.ticks.y  = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    legend.title = element_blank(),
    legend.key= element_rect(colour = "black"),
    aspect.ratio=1) 
  
  
  # Plot V ----------------------------------------------------------------------
  
  # Gather but keep factor order for the columns of V
  df.V <- tidyr::gather(df.V, key='collab', value='value' , 
                        -rowlab, factor_key = TRUE)
  
  df.V$value <- as.factor(df.V$value)
  p <- ggplot2::ggplot(df.V, aes(x=collab, y=rowlab)) + 
    geom_tile(aes(fill = value), colour = "white") + 
    scale_fill_manual(values=c("green", "blue", "red"), na.value = 'white',
                      labels = c("Sí", "Abstención", "No", "Ausente")) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    theme_bw() +  xlab("votaciones") + ylab("") + mytheme +   coord_fixed()
  print(p)
  
  p
}


plot_V <- function(V, rowlabels = FALSE, collabels=FALSE,
                      xlab = "", ylab = "", aspect.ratio = 1){

  base_size <- 7
  axis.text.y <- ggplot2::element_text(angle = 0, size=base_size*1, hjust = 1, 
                                       colour = "black")
  axis.text.x <- ggplot2::element_text(angle = 45, size=base_size*1, hjust = 1, 
                                       colour = "black")
  
  if(is.null(rownames(V))){
    rownames(V) = 1:dim(V)[1]
  }
  
  if (rowlabels == FALSE){ 
    axis.text.y  <- ggplot2::element_blank()
    axis.ticks.y <- ggplot2::element_blank()
  }
  
  if (collabels == FALSE){ 
    axis.text.x  <- ggplot2::element_blank() 
    axis.ticks.x <- ggplot2::element_blank()
  }
  
  # Plot V ----------------------------------------------------------------------
  df.V <- data.frame(V)

  
  #colnames(df.V) <- 1:dim(V)[2]
  # Labs as factors preserving original order
  df.V$rowlab <- rownames(V)
  df.V$rowlab <- factor(df.V$rowlab, levels = df.V$rowlab)
  
  mytheme <- ggplot2::theme(
    axis.text.x = axis.text.x,
    axis.text.y = axis.text.y,
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    aspect.ratio=aspect.ratio)

  
  # Gather but keep factor order for the columns of V
  df.V <- tidyr::gather(df.V, key='collab', value='value' , 
                        -rowlab, factor_key = TRUE)
  

  p <- ggplot2::ggplot(df.V, aes(x=collab, y=rowlab)) + 
    geom_tile(aes(fill = value)) + 
    scale_fill_gradient(low = "white", high = "black", na.value = 'red') +
    #geom_tile(aes(fill = as.factor(value))) + 
    #scale_fill_manual(values=c("white", "black", "blue"), na.value = 'red') +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    theme_bw() +  xlab(xlab) + ylab(ylab) + mytheme +   coord_fixed() 
  print(p)
  
  p
}

plot_dictionary <- function(W, Kmax=100, rowlabels=FALSE, 
                            rowlab_by = 1, base_size = 5,
                            sort=TRUE, rotated=FALSE, aspect.ratio = 1, ...){

  if(hasArg(labels)){
    labels <- list(...)$labels
    rownames(W) <- labels
  }
  
  
  if(sort){
    idx.sorted.W <- order(-colSums(W)) # columns sort by size in W
    W <- W[,idx.sorted.W]
  } 
  
  df.W <- data.frame(W)
  colnames(df.W) <- 1:dim(W)[2]
  
  if(is.null(rownames(W)) && rowlabels){ 
    warning("W has no rownames. Using labels = FALSE")
    rowlabels=FALSE
  }
  
  if(rowlabels){
    df.W$rowlab <- rownames(W)
    mytheme <- ggplot2::theme(
      panel.grid.minor = element_blank(),
      axis.text.y = ggplot2::element_text(size=base_size, hjust = 1, 
                           colour = "black"),
      legend.position = "none",
      plot.margin=unit(c(0, 0,0,-3), "mm"),
      aspect.ratio=aspect.ratio)
  }
  if(!rowlabels){
    df.W$rowlab <- 1:dim(W)[1]
    mytheme <- ggplot2::theme(
                      axis.text.y = ggplot2::element_blank(),
                      axis.ticks.y = ggplot2::element_blank(),
                      panel.grid.minor = ggplot2::element_blank(),
                      legend.position = "none",
                      plot.margin=unit(c(0, 0,0,-5), "mm"),
                      aspect.ratio=aspect.ratio)
  }
  df.W$rowlab <- factor(df.W$rowlab, levels = df.W$rowlab)
  df.W <- tidyr::gather(df.W, key='dimension', value='value' , -rowlab)
  df.W$dimension <- as.numeric(df.W$dimension)
  
  rowbreaks <- levels(df.W$rowlab)[seq(1,length(levels(df.W$rowlab)), 
                                         by=rowlab_by)]

  
  if(!rotated){
    p <- ggplot2::ggplot(df.W %>% filter(dimension <Kmax), 
                         aes(x=factor(dimension), y=rowlab)) + 
      geom_tile(aes(fill = value), colour = "white") + 
      scale_fill_gradient(low = "white", high = "black", na.value = 'red') +
      scale_x_discrete(expand = c(0, 0)) +
      scale_y_discrete(expand = c(0, 0), breaks = rowbreaks) +
      xlab("component") +
      ylab("") +
      theme_bw() + mytheme +   coord_fixed() 
    print(p)
  } else{
    # Features in horizontal
    mytheme <- ggplot2::theme(
      panel.grid.minor = element_blank(),
      legend.position = "none",
      axis.text.x = ggplot2::element_text(angle = 90, size=base_size*1, hjust = 1, 
                            colour = "black"),
      axis.text.y = ggplot2::element_text(size=base_size*1, hjust = 1, 
                                           colour = "black"),
      plot.margin=unit(c(0,1,0,-3), "mm"),
      aspect.ratio=1/15)

    # if UN dataset, use the country codes
    #df$country_code[match(rownames(E_W), un_votes$country)]
    #df.W$rowlab <- un_votes$country_code[match(df.W$rowlab, un_votes$country)]
    
    p <- ggplot2::ggplot(df.W %>% filter(dimension <Kmax), 
                         aes(y=factor(dimension), x=rowlab)) + 
      geom_tile(aes(fill = value), colour = "white") + 
      scale_fill_gradient(low = "white", high = "black", na.value = 'red') +
      scale_x_discrete(expand = c(0, 0)) +
      scale_y_discrete(expand = c(0, 0)) +
      xlab("") +
      ylab("") +
      theme_bw() + mytheme
    print(p)
  }
  p
}

plot_multiple_dictionaries <- function(list.W, names , Kmax=100, rowlabels=FALSE, 
                              sort=TRUE, rotated=FALSE, aspect.ratio = 1, 
                              idx.subset = NULL,
                              ...){
  
  # Convert dictionaries to dataframe
  if(length(list.W) != length(names)) stop("names do not match dictionaries")
  
  base_size <- 10
  
  if(hasArg(labels)){
    labels <- list(...)$labels
    for(i in 1:length(list.W)){
      rownames(list.W[[i]]) <- labels
    }
    
  }
  
  if(!is.null(idx.subset)){
    #idx.active <- sample(dim(list.W[[1]])[1])[1:100]  # random
    for(i in 1:length(list.W)){
      list.W[[i]] <- list.W[[i]][idx.subset,]
    }
  }
  
  
  # Sort according to the first dictionary
  if(sort){
    for(i in 1:length(list.W)){
      idx.sorted.W <- order(-colSums(list.W[[i]]))
      list.W[[i]] <- list.W[[i]][,idx.sorted.W]
    }
  } 
  
  rownames <- rownames(list.W[[1]])
  df.W <- data.frame()
  for(i in 1:length(list.W)){
    W <- list.W[[i]]
    if(is.null(rownames(W))) {
      rownames(W) <- 1:nrow(W)
    }
    df <- as.data.frame(W)
    names(df) <- 1:ncol(df)
    df$rowname <- rownames(df)
    df.W_ <- tidyr::gather(df, key='dimension', value='value', -rowname)
    df.W_$name <- names[i]
    df.W <- bind_rows(df.W, df.W_)
  }
  df.W$rowname <- factor(df.W$rowname, levels = rownames)
  df.W$name <- factor(df.W$name, levels=names)
  
  if(rowlabels){
    mytheme <- ggplot2::theme(
      strip.background =element_rect(fill="white"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      legend.position = "none",
      plot.margin=unit(c(0, 0,0,-3), "mm"),
      #plot.background = element_rect(fill="red"),
      axis.text.y = ggplot2::element_text(size=base_size-1, hjust = 1, 
                                          colour = "black"),
      panel.background=element_rect(fill="lightgrey", colour="red", linetype = "blank"),
      aspect.ratio=aspect.ratio)
  }
  if(!rowlabels){
    mytheme <- ggplot2::theme(
      strip.background =element_rect(fill="white"),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "none",
      plot.margin=unit(c(0, 0,0,-5), "mm"),
      #plot.background = element_rect(fill="red"),
      aspect.ratio=aspect.ratio)
  }
  df.W$dimension <- as.numeric(df.W$dimension)
  
  if(!rotated){
    p <- ggplot2::ggplot(df.W %>% filter(dimension <= Kmax), 
                         aes(x=factor(dimension), y=rowname)) + 
      geom_tile(aes(fill = value), colour = "white") + 
      scale_fill_gradient(low = "white", high = "black", na.value = 'red') +
      scale_x_discrete(expand = c(0, 0)) +
      scale_y_discrete(expand = c(0, 0)) +
      facet_grid(. ~ name) +
      xlab("component") +
      ylab("") +
      theme_bw() + mytheme +   coord_fixed() 
    print(p)
  } else{
    # Features in horizontal
    mytheme <- ggplot2::theme(
      panel.grid.minor = element_blank(),
      strip.background =element_rect(fill="white"),
      legend.position = "none",
      axis.text.x = ggplot2::element_text(angle = 90, size=base_size*1, hjust = 1, 
                                          colour = "black"),
      axis.text.y = ggplot2::element_text(size=base_size*1, hjust = 1, 
                                          colour = "black"),
      plot.margin=unit(c(0,1,0,-3), "mm"),
      aspect.ratio=aspect.ratio)
    
    # if UN dataset, use the country codes
    #df$country_code[match(rownames(E_W), un_votes$country)]
    #df.W$rowlab <- un_votes$country_code[match(df.W$rowlab, un_votes$country)]
    
    p <- ggplot2::ggplot(df.W %>% filter(dimension <= Kmax), 
                         aes(y=factor(dimension), x=rowname)) + 
      geom_tile(aes(fill = value), colour = "white") + 
      scale_fill_gradient(low = "white", high = "black", na.value = 'red') +
      scale_x_discrete(expand = c(0, 0)) +
      scale_y_discrete(expand = c(0, 0)) +
      xlab("") +
      ylab("") +
      facet_grid(name ~ .) +
      theme_bw() + mytheme #+   coord_fixed() 
    print(p)
  }
  p
}


#' @title Plot trace of W dimensions
plot_trace_W_dimensions <- function(W_samples){
  J <- dim(W_samples)[3]
  K <- dim(W_samples)[2]
  
  dimensions <- t(apply(W_samples, 3, function(W) {colSums(W)}))
  df.dimensions <- as.data.frame(dimensions)
  names(df.dimensions) <- 1:K
  df.dimensions$sample <- 1:J
  df.dimensions <- gather(df.dimensions, key = "k", value = "norm", -sample)

  p <- ggplot(df.dimensions, aes(x=sample, y =norm, group = k, color= k))+
    geom_line() + 
    theme_bw() +
    theme(legend.position = "none") +
    ylab("k")
  print(p)
  p
}


#' #' @title Plot trace of w_k * h_k
#' plot_trace_WH_dimensions <- function(W_samples, H_samples){
#'   J <- dim(W_traces)[3]
#'   K <- dim(W_samples)[2]
#'   m <- matrix(NA, J, K)
#'   for(j in 1:J){
#'     W_samples[,,j] %*% H_samples[,,j]
#'   }
#'   F <- dim(Z_samples)[1]
#'   N <- dim(Z_samples)[2]
#'   
#'   
#'   df.traces <- data.frame(sample = 1:J)
#'   for(i in 1:nelements){
#'     df.traces <- cbind(df.traces, Z_samples[f_selected[i], n_selected[i],])
#'   }
#'   names(df.traces)[2:(nelements+1)] <- 1:nelements
#'   df.traces <- gather(df.traces, key="element", value = "k", -sample)
#'   p <- ggplot(df.traces, aes(x=sample, y = k, group = element, color= element))+
#'     facet_grid(element~.)+ 
#'     geom_line() + 
#'     theme_bw() +
#'     theme(legend.position = "none")+
#'     ylab("k")
#'   p
#' }


#' @title Plot number of active components during the Gibbs sampler
plot_trace_size_components <- function(model){
  Z_samples <- model$Z_samples
  J <- dim(Z_samples)[3]
  Kmax <- max(Z_samples, na.rm = TRUE)
  df.components <- as.data.frame(t(apply(Z_samples, 3, function(Z) {
                            tabulate(c(Z)+1, nbins=Kmax)
                            })))
  names(df.components) <- 1:Kmax
  df.components$sample <- 1:J
  df <- gather(df.components, key = "component", value = "size", -sample)
  ggplot(df, aes(x=sample, y = size, group = component)) + 
    geom_line() + 
    theme_bw() + 
    ylab("components size")
}

#' @title Plot trace of some elements Z(f,n)
plot_trace_elements <- function(model, nelements = 10){
  Z_samples <- model$Z_samples
  F <- dim(Z_samples)[1]
  N <- dim(Z_samples)[2]
  J <- dim(Z_samples)[3]
  f_selected <- sample(F, nelements, replace=FALSE)
  n_selected <- sample(N, nelements, replace=FALSE)
  Z_samples[1,1,]
  
  df.traces <- data.frame(sample = 1:J)
  for(i in 1:nelements){
    df.traces <- cbind(df.traces, Z_samples[f_selected[i], n_selected[i],])
  }
  names(df.traces)[2:(nelements+1)] <- 1:nelements
  df.traces <- gather(df.traces, key="element", value = "k", -sample)
  p <- ggplot(df.traces, aes(x=sample, y = k, group = element, color= element))+
      facet_grid(element~.)+ 
      geom_line() + 
      theme_bw() +
      theme(legend.position = "none")+
      ylab("k")
  p
}

#' @title Plot number of active components during the Gibbs sampler
plot_trace_ncomponents <- function(model){
  Z_samples <- model$Z_samples
  J <- dim(Z_samples)[3]
  ncomponents <- apply(Z_samples, 3, function(Z) {length(unique(c(Z)))})
  df <- data.frame(sample = 1:J, k=ncomponents)
  ggplot(df, aes(x=sample, y = k)) + geom_point() + 
    theme_bw() + 
    ylab("active components")
}


#' @title Plot trace of the Variational Lower Bound
plot_trace_lowerbound <- function(model){
  lb <- c(model$lowerbounds)
  df <- data.frame(iter = 1:length(lb), lb=lb)
  p <- ggplot(df, aes(x=iter, y=lb)) + geom_line(size=0.5) + 
    #geom_point(size=0.2)+
    theme_bw() +
    theme(axis.text = element_text(color = "black")) +
    ylab("Lower Bound") + xlab("iteration")
  print(p)
  p
}