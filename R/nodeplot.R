### function node_plot plots the data assigned to a specified node, projected into a
### two dimsional subspace, to illustrate the split at that node (or the would-be split for a leaf node)
## arguments:
# sol = cluster solution arising from any of the clustering algorithms in the package
# node = either the node number to be viewed (nodes are listed in the order
#        they are added to the model), or a vector of length 2 specifying the
#        depth and position of the node within the hierarchy.
# labels = vector of length n. If class labels are given then the
#          performance at the specified node is given.

node_plot <- function(sol, node, labels = NULL){
  op <- par(no.readonly = TRUE)

  # if node location is given by its position (rather than number), then determine its number

  if(length(node)>1){
    for(i in 1:length(sol$Nodes)){
      if(sum(node==sol$Nodes[[i]]$location)==2){
        node = i
        break
      }
    }
    if(length(node)>1) stop('You must specify a node location within the given hierarchy')
  }

  X <- sol$data

  # if labels are given, then ensure they are integer valued for the purpose of colour plots

  if(!is.null(labels)){
    lab_new <- numeric(length(labels))
    u <- unique(labels)
    for(i in 1:length(u)) lab_new[which(labels==u[i])] = i
    labels <- lab_new
  }

  # determine dimensions and layout for the plot

  d <- max(sol$model[,1])
  w <- subtree_width(sol, 1)
  l.mat <- matrix(0, 4, 3)
  l.mat[1,] <- c(2, 2, 1)
  l.mat[2,] <- c(2, 2, 1)
  l.mat[3,] <- c(2, 2, 3)
  l.mat[4,] <- c(2, 2, 3)
  layout(l.mat)

  # plot the full hierarchical structure in the upper corner

  par(mar = c(0, 0, 0, 0))
  plot(0, cex = 0, ylim = c(0, 1), xlim = c(0, 1), xaxt = 'n', yaxt = 'n', bty = 'n')
  for(i in 1:d){
    ixs = which(sol$model[,1]==i)
    for(j in ixs){
      loc <- sol$model[j,]
      is.leaf = 1-sum((sol$model[,1]==(i+1))*(sol$model[,2]==(2*loc[2])))
      if(is.leaf){
        y1 <- 1-(i-1)/d
        y2 <- 1-i/d
        x1 <- (2*loc[2]-1)/2^loc[1]
        x2 <- (2*loc[2]-1)/2^loc[1]
        segments(x1, y1, x2, y2)
      }
      else{
        y1 <- 1-(i-1)/d
        y2 <- 1-i/d
        x1 <- (2*loc[2]-1)/2^loc[1]
        x2 <- (2*loc[2]-1)/2^loc[1]
        segments(x1, y1, x2, y2)
        x1 <- (2*(2*loc[2]-1)-1)/2^(loc[1]+1)
        x2 <- (2*(2*loc[2])-1)/2^(loc[1]+1)
        segments(x1, y2, x2, y2)
      }
      if(j==node) points((2*loc[2]-1)/2^loc[1], 1-i/d, col = rgb(1, 0, 0, .5), pch = 16, cex = 3)
    }
  }


  # plot the projected data in the main body of the plot

  par(mar = c(2, 2, 0, 2))

  is.leaf = 1-sum((sol$model[,1]==(sol$model[node,1]+1))*(sol$model[,2]==(2*sol$model[node,2])))

  if(ncol(X)>2) v2 <- rARPACK::eigs_sym(cov(X[sol$Nodes[[node]]$ixs,]-X[sol$Nodes[[node]]$ixs,]%*%sol$Nodes[[node]]$v%*%t(sol$Nodes[[node]]$v)), 1)$vectors[,1]
  else v2 <- eigen(cov(X[sol$Nodes[[node]]$ixs,]-X[sol$Nodes[[node]]$ixs,]%*%sol$Nodes[[node]]$v%*%t(sol$Nodes[[node]]$v)))$vectors[,1]

  Xp <- X[sol$Nodes[[node]]$ixs,]%*%cbind(sol$Nodes[[node]]$v, v2)
  if(sol$method=='MDH') den <- density(Xp[,1], bw = sol$Nodes[[node]]$params$h)
  else den <- density(Xp[,1])

  if(is.null(labels)){
    if(is.leaf){
      if(sol$model[node,2]%%2) plot(Xp, col = 4, tck = .02, yaxt = 'n', bty = 'n')
      else plot(Xp, col = 2, tck = .02, yaxt = 'n', bty = 'n')
    }
    else plot(Xp, col = (Xp[,1]<sol$Nodes[[node]]$b)*2+2, tck = .02, yaxt = 'n', bty = 'n')
  }
  else{
    plot(Xp, col = labels[sol$Nodes[[node]]$ixs], tck = .02, yaxt = 'n', bty = 'n')
  }
  axis(2, labels = round(seq(min(Xp[,2]), min(Xp[,2])+(max(Xp[,2])-min(Xp[,2]))*450/512, length = 7), 1), at = seq(min(Xp[,2]), min(Xp[,2])+(max(Xp[,2])-min(Xp[,2]))*450/512, length = 7), tck = .02)
  lines(den$x, den$y/max(den$y)*(max(Xp[,2])-min(Xp[,2]))+min(Xp[,2]), lwd = 2)
  axis(4, labels = round(seq(0, max(den$y)*450/512, length = 7), 2), at = seq(min(Xp[,2]), min(Xp[,2])+(max(Xp[,2])-min(Xp[,2]))*450/512, length = 7), tck = .02)
  abline(h = min(Xp[,2]))

  if(!is.null(labels)){
    T = table(Xp[,1]<sol$Nodes[[node]]$b, labels[sol$Nodes[[node]]$ixs])
    split = apply(T, 2, which.max)

    if(max(split)==1){
      lines(den$x, den$y/max(den$y)*(max(Xp[,2])-min(Xp[,2]))+min(Xp[,2]), col = 2, lwd = 2)
    }
    else if(min(split)==2){
      lines(den$x, den$y/max(den$y)*(max(Xp[,2])-min(Xp[,2]))+min(Xp[,2]), col = 4, lwd = 2)
    }
    else{
      ixl <- which(labels[sol$Nodes[[node]]$ixs]%in%sort(unique(labels[sol$Nodes[[node]]$ixs]))[which(split==2)])
      ixr <- which(labels[sol$Nodes[[node]]$ixs]%in%sort(unique(labels[sol$Nodes[[node]]$ixs]))[which(split==1)])
      denl <- density(Xp[ixl,1], bw = den$bw)
      denr <- density(Xp[ixr,1], bw = den$bw)
      lines(denl$x, denl$y*length(ixl)/length(sol$Nodes[[node]]$ixs)/max(den$y)*(max(Xp[,2])-min(Xp[,2]))+min(Xp[,2]), col = 4, lwd = 2)
      lines(denr$x, denr$y*length(ixr)/length(sol$Nodes[[node]]$ixs)/max(den$y)*(max(Xp[,2])-min(Xp[,2]))+min(Xp[,2]), col = 2, lwd = 2)
    }
  }


  # add details of the split

  if(is.leaf){
    abline(v = sol$Nodes[[node]]$b, col = rgb(.5, .5, .5), lty = 2, lwd = 2)
    par(mar = c(0, 0, 0, 0))
    plot(.05, cex = 0, ylim = c(0, 1), xlim = c(0, 1), bty = 'n', xaxt = 'n', yaxt = 'n')
    text(.05, .95, paste('node: ', as.character(node)), pos = 4, offset = 0)
    text(.05, .85, paste('depth: ', as.character(sol$model[node,1])), pos = 4, offset = 0)
    text(.05, .75, paste('position: ', as.character(sol$Nodes[[node]]$location[2])), pos = 4, offset = 0)
    text(.05, .65, paste('n: ', as.character(length(sol$Nodes[[node]]$ixs))), pos = 4, offset = 0)
    if(sol$method=='MCDC') text(.05, .55, paste('variance ratio: ', substr(as.character(sol$Nodes[[node]]$fval), 1, 5)), pos = 4, offset = 0)
    else if(sol$method=='MDH') text(.05, .55, paste('relative depth: ', substr(as.character(sol$Nodes[[node]]$rel.dep), 1, 5)), pos = 4, offset = 0)
    else if(sol$method=='NCutH') text(.05, .55, paste('normalised cut: ', substr(as.character(sol$Nodes[[node]]$fval), 1, 5)), pos = 4, offset = 0)
    if(!is.null(labels)){
      prf = max(table(labels[sol$Nodes[[node]]$ixs]))/length(sol$Nodes[[node]]$ixs)
      text(.05, .45, paste('cluster purity: ', substr(as.character(prf), 1, 5)), pos = 4, offset = 0)
    }
    else text(.05, .45, 'cluster purity: -', pos = 4, offset = 0)
  }
  else{
    abline(v = sol$Nodes[[node]]$b, col = 2, lwd = 2)
    par(mar = c(0, 0, 0, 0))
    plot(0, cex = 0, ylim = c(0, 1), xlim = c(0, 1), bty = 'n', xaxt = 'n', yaxt = 'n')
    text(.05, .95, paste('node: ', as.character(node)), pos = 4, offset = 0)
    text(.05, .85, paste('depth: ', as.character(sol$model[node,1])), pos = 4, offset = 0)
    text(.05, .75, paste('position: ', as.character(sol$Nodes[[node]]$location[2])), pos = 4, offset = 0)
    text(.05, .65, paste('n: ', as.character(length(sol$Nodes[[node]]$ixs))), pos = 4, offset = 0)
    if(sol$method=='MCDC') text(.05, .55, paste('variance ratio: ', substr(as.character(sol$Nodes[[node]]$fval), 1, 5)), pos = 4, offset = 0)
    else if(sol$method=='MDH') text(.05, .55, paste('relative depth: ', substr(as.character(sol$Nodes[[node]]$rel.dep), 1, 5)), pos = 4, offset = 0)
    else if(sol$method=='NCutH') text(.05, .55, paste('normalised cut: ', substr(as.character(sol$Nodes[[node]]$fval), 1, 5)), pos = 4, offset = 0)
    if(!is.null(labels)){
      prf = success_ratio((X[sol$Nodes[[node]]$ixs,]%*%sol$Nodes[[node]]$v<sol$Nodes[[node]]$b), labels[sol$Nodes[[node]]$ixs])
      text(.05, .45, paste('success ratio: ', substr(as.character(prf[1]), 1, 5)), pos = 4, offset = 0)
    }
    else text(.05, .45, 'success ratio: -', pos = 4, offset = 0)
  }
  par(op)
}

### function hp_plot provides an illustration of a single hyperplane partition, arising from
### one of mdh, mch and ncuth
## arguments:
# sol = solution from one of mdh, mch and ncuth
# X = data matrix used to generate the partition
# labels = vector of class labels. If provided then points are plotted with their class label's colour

## mch, mdh and ncuth provide a list of potential splits, arising from each of the initialisatinos provided
## the default is to plot the solution with the optimal value of its projection index. If one wishes to
## specify which hyperplane to visualise, use hp_plot(sol[[i]], X) for the i-th solution.

hp_plot <- function(sol, X, labels = NULL){

  op <- par(no.readonly = TRUE)

  par(mar = c(2, 2, 2, 2))

  # if labels are given, then ensure they are integer valued for the purpose of colour plots

  if(!is.null(labels)){
    lab_new <- numeric(length(labels))
    u <- unique(labels)
    for(i in 1:length(u)) lab_new[which(labels==u[i])] = i
    labels <- lab_new
  }

  # if multiple hyperplanes are present in the solution then select the best

  if(is.null(sol$fval)){
    if(sol[[1]]$method=='MCDC') ix.opt <- which.max(unlist(lapply(sol, function(hyp) hyp$fval)))
    else ix.opt <- which.min(unlist(lapply(sol, function(hyp) hyp$fval)))
    sol <- sol[[ix.opt]]
  }

  # project data into two-dimensional subspace for plotting
  if(ncol(X)>2) v2 <- rARPACK::eigs_sym(cov(X-X%*%sol$v%*%t(sol$v)), 1)$vectors[,1]
  else v2 <- eigen(cov(X-X%*%sol$v%*%t(sol$v)))$vectors[,1]

  Xp <- X%*%cbind(sol$v, v2)

  # compute the external quality of the split through success ratio (if labels are provided)

  if(!is.null(labels)){
    prf <- success_ratio((Xp[,1]<sol$b), labels)
  }
  else prf <- '-'


  # plot the projected points and the estimated density along the optimal projection

  if(sol$method=='MDH') den <- density(Xp[,1], bw = sol$params$h)
  else den <- density(Xp[,1])

  if(is.null(labels)){
    plot(Xp, col = (Xp[,1]<sol$b)*2+2, tck = .02, yaxt = 'n')
  }
  else{
    plot(Xp, col = labels, tck = .02, yaxt = 'n')
  }
  abline(v = sol$b, col = 2, lwd = 2)
  axis(2, labels = round(seq(min(Xp[,2]), min(Xp[,2])+(max(Xp[,2])-min(Xp[,2]))*450/512, length = 7), 1), at = seq(min(Xp[,2]), min(Xp[,2])+(max(Xp[,2])-min(Xp[,2]))*450/512, length = 7), tck = .02)
  lines(den$x, den$y/max(den$y)*(max(Xp[,2])-min(Xp[,2]))+min(Xp[,2]), lwd = 2)
  axis(4, labels = round(seq(0, max(den$y)*450/512, length = 7), 2), at = seq(min(Xp[,2]), min(Xp[,2])+(max(Xp[,2])-min(Xp[,2]))*450/512, length = 7), tck = .02)
  abline(h = min(Xp[,2]))

  if(!is.null(labels)){
    T = table(Xp[,1]<sol$b, labels)
    split = apply(T, 2, which.max)

    if(max(split)==1){
      lines(den$x, den$y/max(den$y)*(max(Xp[,2])-min(Xp[,2]))+min(Xp[,2]), col = 2, lwd = 2)
    }
    else if(min(split)==2){
      lines(den$x, den$y/max(den$y)*(max(Xp[,2])-min(Xp[,2]))+min(Xp[,2]), col = 4, lwd = 2)
    }
    else{
      ixl <- which(labels%in%sort(unique(labels))[which(split==2)])
      ixr <- which(labels%in%sort(unique(labels))[which(split==1)])
      denl <- density(Xp[ixl,1], bw = den$bw)
      denr <- density(Xp[ixr,1], bw = den$bw)
      lines(denl$x, denl$y*length(ixl)/length(labels)/max(den$y)*(max(Xp[,2])-min(Xp[,2]))+min(Xp[,2]), col = 4, lwd = 2)
      lines(denr$x, denr$y*length(ixr)/length(labels)/max(den$y)*(max(Xp[,2])-min(Xp[,2]))+min(Xp[,2]), col = 2, lwd = 2)
    }
  }

  # add details of the solutions

  if(sol$method=='MCDC') title(main = paste('variance ratio: ', substr(as.character(sol$fval), 1, 5), 'success ratio: ', substr(as.character(prf[1]), 1, 5)))
  else if(sol$method=='MDH') title(main = paste('relative depth: ', substr(as.character(sol$rel.dep), 1, 5), 'success ratio: ', substr(as.character(prf[1]), 1, 5)))
  else if(sol$method=='NCutH') title(main = paste('normalised cut: ', substr(as.character(sol$fval), 1, 5), 'success ratio: ', substr(as.character(prf[1]), 1, 5)))

  par(op)
}
