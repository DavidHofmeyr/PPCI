### plot functions for classes 'ppci_cluster_solution' and 'ppci_hyperplane_solution' provide
### wrappers for functions tree_plot, node_plot and hp_plot. Plotting an object of class
### 'ppci_cluster_solution' will use tree_plot, representing the entire solution, if the node
### number is not specified, otherwise it will use the function node_plot which provides a more
### detailed view of the chosen node in the hierarchical model. plot function for class
### 'ppci_projection_solution' plots the first two dimensions of the projected data. If argument
### pairs is given, then pairs of dimensions are plotted agains one another.

plot.ppci_cluster_solution <- function(x, node = NULL, labels = NULL, node.numbers = NULL, transparency = NULL, ...){
  #control <- list(...)
  if(is.null(node)) tree_plot(x, labels, node.numbers, transparency)
  else node_plot(x, node, labels, transparency)
}

plot.ppci_hyperplane_solution <- function(x, labels = NULL, transparency = NULL, ...){
  #control <- list(...)
  hp_plot(x, labels, transparency)
}

plot.ppci_projection_solution <- function(x, labels = NULL, pairs = NULL, PCA = NULL, transparency = NULL, ...){
  #control <- list(...)
  if(is.null(transparency)) transparency = 0
  if(!is.null(PCA) && PCA) x$fitted <- x$fitted%*%eigen(cov(x$fitted))$vectors
  if(is.null(pairs)){
    if(is.null(labels)){
      if(transparency==0) plot(x$fitted[,1:2])
      else plot(x$fitted[,1:2], col = rgb(0, 0, 0, transparency), pch = 16)
    }
    else{
      lab_new <- numeric(length(labels))
      u <- unique(labels)
      for(i in 1:length(u)) lab_new[which(labels==u[i])] = i
      labels <- lab_new
      if(transparency==0) plot(x$fitted[,1:2], col = labels)
      else plot(x$fitted[,1:2], col = sapply(labels, function(c){
        cl = col2rgb(c)
        rgb(cl[1]/255, cl[2]/255, cl[3]/255, transparency)
      }), pch = 16)
    }
  }
  else if(is.numeric(pairs)){
    if(is.null(labels)){
      if(transparency==0) pairs(x$fitted[,1:pairs])
      else pairs(x$fitted[,1:pairs], col = rgb(0, 0, 0, transparency), pch = 16)
    }
    else{
      lab_new <- numeric(length(labels))
      u <- unique(labels)
      for(i in 1:length(u)) lab_new[which(labels==u[i])] = i
      labels <- lab_new
      if(transparency==0) pairs(x$fitted[,1:pairs], col = labels)
      else pairs(x$fitted[,1:pairs], col = sapply(labels, function(c){
        cl = col2rgb(c)
        rgb(cl[1]/255, cl[2]/255, cl[3]/255, transparency)
      }), pch = 16)
    }
  }
  else stop('pairs must be a positive integer')
}

### function tree_plot provides an illustration of a complete hierarchical clustering solution
### arising from any of the clustering algorithms in the package. each partition in the slution
### is visualised through a two-dimensional plot of the data being split at the relevant node.
## arguments:
# sol = clustering solution arising from one of the methods in the package
# labels = vector of class labels. If provided then plots indicate the class membership of the data. Otherwise
#           points are coloured according to the different two-way partitions
# node.numbers = logical. should the order in which nodes were added to the solution be indicated on the plot or not.

tree_plot <- function(sol, labels = NULL, node.numbers = NULL, transparency = NULL){
  if(is.null(transparency)) transparency = 0
  if(is.null(node.numbers)) node.numbers <- TRUE
  op <- par(no.readonly = TRUE)
  par(mar = c(0, 0, 0, 0))
  X <- sol$data
  plot(0, xlim = c(0, 1), ylim = c(0, 1), cex = 0, xaxt = 'n', yaxt = 'n')

  # determine the geometry of the tree for plotting

  d <- max(sol$model[,1])
  ns <- sapply(1:d, function(i) sum(sol$model[,1]==i))
  w <- max(ns)
  width <- 1/w
  height <- .8/d

  # if labels are supplied, ensure they are integers 1:K (K the number of classes)

  if(!is.null(labels)){
    lab_new <- numeric(length(labels))
    u <- unique(labels)
    for(i in 1:length(u)) lab_new[which(labels==u[i])] = i
    labels <- lab_new
  }

  # add the plots of each partition in the hierarchy recursively using add_subtree

  add_subtree(sol, X, 1, 0, 1, 1, height, width, labels, node.numbers, transparency)
  par(op)
}

### function add_subtree recursively adds plots of individual two-way partitions to the existing plot.
### function is not intended to be called outside of tree_plot

add_subtree <- function(sol, X, node, L, U, y, h, w, labels, node.numbers, transparency){
  is.leaf <- 1-sum((sol$model[,1]==(sol$model[node,1]+1))*(sol$model[,2]==sol$model[node,2]*2))
  if(is.leaf){
    if(ncol(X)>2) v2 <- rARPACK::eigs_sym(cov(X[sol$Nodes[[node]]$ixs,]-X[sol$Nodes[[node]]$ixs,]%*%sol$Nodes[[node]]$v%*%t(sol$Nodes[[node]]$v)), 1)$vectors[,1]
    else v2 <- eigen(cov(X[sol$Nodes[[node]]$ixs,]-X[sol$Nodes[[node]]$ixs,]%*%sol$Nodes[[node]]$v%*%t(sol$Nodes[[node]]$v)))$vectors[,1]
    Xp <- X[sol$Nodes[[node]]$ixs,]%*%cbind(sol$Nodes[[node]]$v, v2)
    Xp[,1] <- w*.7*(Xp[,1]-min(Xp[,1]))/(max(Xp[,1])-min(Xp[,1])) + (L+U-w*.7)/2
    Xp[,2] <- h*(Xp[,2]-min(Xp[,2]))/(max(Xp[,2])-min(Xp[,2])) + y - h
    if(is.null(labels)){
      if(sol$model[node,2]%%2){
        if(transparency==0) points(Xp, cex = .5, col = 4)
        else points(Xp, cex = .5, col = rgb(0, 0, 1, transparency), pch = 16)
      }
      else{
        if(transparency==0) points(Xp, col = 2, cex = .5)
        else points(Xp, col = rgb(1, 0, 0, transparency), pch = 16, cex = .5)
      }
    }
    else{
      if(transparency==0) points(Xp, col = labels[sol$Nodes[[node]]$ixs], cex = .5)
      else points(Xp, col = sapply(labels[sol$Nodes[[node]]$ixs], function(c){
        cl = col2rgb(c)
        rgb(cl[1]/255, cl[2]/255, cl[3]/255, transparency)
      }), pch = 16, cex = .5)
    }

    if(node.numbers) text((L+U-w*.8)/2, y-h*.9, as.character(node), cex = .7)
  }
  else{
    if(ncol(X)>2) v2 <- rARPACK::eigs_sym(cov(X[sol$Nodes[[node]]$ixs,]-X[sol$Nodes[[node]]$ixs,]%*%sol$Nodes[[node]]$v%*%t(sol$Nodes[[node]]$v)), 1)$vectors[,1]
    else v2 <- eigen(cov(X[sol$Nodes[[node]]$ixs,]-X[sol$Nodes[[node]]$ixs,]%*%sol$Nodes[[node]]$v%*%t(sol$Nodes[[node]]$v)))$vectors[,1]
    Xp <- X[sol$Nodes[[node]]$ixs,]%*%cbind(sol$Nodes[[node]]$v, v2)
    Xp[,1] <- w*(Xp[,1]-min(Xp[,1]))/(max(Xp[,1])-min(Xp[,1])) + (L+U-w)/2
    Xp[,2] <- h*(Xp[,2]-min(Xp[,2]))/(max(Xp[,2])-min(Xp[,2])) + y - h
    if(is.null(labels)){
      cols <- numeric(nrow(Xp)) + 4
      cols[which(X[sol$Nodes[[node]]$ixs,]%*%sol$Nodes[[node]]$v>sol$Nodes[[node]]$b)] <- 2
      if(transparency==0) points(Xp, col = cols, cex = .5)
      else points(Xp, col = sapply(cols, function(c){
        cl = col2rgb(c)
        rgb(cl[1]/255, cl[2]/255, cl[3]/255, transparency)
      }), pch = 16, cex = .5)
    }
    else{
      if(transparency==0) points(Xp, col = labels[sol$Nodes[[node]]$ixs], cex = .5)
      else points(Xp, col = sapply(labels[sol$Nodes[[node]]$ixs], function(c){
        cl = col2rgb(c)
        rgb(cl[1]/255, cl[2]/255, cl[3]/255, transparency)
      }), pch = 16, cex = .5)
    }

    if(node.numbers) text((L+U-w*1.1)/2, y-h*.9, as.character(node), cex = .7)

    k1 <- which(((sol$model[,1]==(sol$model[node,1]+1))*(sol$model[,2]==(sol$model[node,2]*2-1)))==1)
    k2 <- which(((sol$model[,1]==(sol$model[node,1]+1))*(sol$model[,2]==sol$model[node,2]*2))==1)
    w1 <- subtree_width(sol, k1)
    w2 <- subtree_width(sol, k2)
    if(w1>w2) w1 <- w1*1.2
    else if(w1<w2) w2 <- w2*1.2
    M <- L + (U-L)*w1/(w1+w2)
    segments((L+U)/2, y-h, (L+U)/2, y-1.125*h)
    segments((L+M)/2, y-1.125*h, (M+U)/2, y-1.125*h)
    segments((L+M)/2, y-1.125*h, (L+M)/2, y-1.25*h)
    segments((U+M)/2, y-1.125*h, (U+M)/2, y-1.25*h)
    add_subtree(sol, X, k1, L, M, y - 1.25*h, h, (M-L)/w1, labels, node.numbers, transparency)
    add_subtree(sol, X, k2, M, U, y - 1.25*h, h, (U-M)/w2, labels, node.numbers, transparency)
  }
}

### function subtree_width used to determine geometry of the plots in plot_tree. Not intended to be
### called directly.

subtree_width <- function(sol, node){
  upper <- lower <- sol$model[node,2]
  depth <- sol$model[node,1]
  rows <- c(node)
  for(i in (depth+1):max(sol$model[,1])){
    upper <- 2*upper
    lower <- 2*lower - 1
    add.rows <- (sol$model[,1]==i)*(sol$model[,2]>=lower)*(sol$model[,2]<=upper)
    rows <- c(rows, which(add.rows==1))
  }
  if(length(rows)==1) return(1)
  sub.model <- sol$model[rows,]
  ns <- sapply(unique(sub.model[,1]), function(i) sum(sub.model[,1]==i))
  max(ns)
}

