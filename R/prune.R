### function tree_prune allows a user to modify an existing hierarchical clustering solution by
### pruning a subtree arising from a specified node.
## arguments:
# sol = clustering solution arising from one of the clustering algorithms in the package
# node = either the node number (based on the order of addition to the model), or a vector specifying the
#       location of the node in the tree (c(depth, position at depth))

## output is a clustering solution after modification. This solution has the same form and type as the original solution

tree_prune <- function(sol, node){

  # if node is specified by location, determine its node number

  if(length(node)>1){
    for(i in 1:length(sol$Nodes)){
      if(sum(sol$Nodes[[i]]$location==node)==2){
        node <- i
        break
      }
    }
    if(length(node)>1) stop('You must specify a node location within the given hierarchy')
  }

  is.leaf = 1-sum((sol$model[,1]==(sol$model[node,1]+1))*(sol$model[,2]==(2*sol$model[node,2])))
  if(is.leaf) stop('pruning from a leaf node will have no effect')

  # find the nodes in the subtree arising from the specified node

  upper <- lower <- sol$model[node,2]
  depth <- sol$model[node,1]
  rows <- c()
  for(i in (depth+1):max(sol$model[,1])){
    upper <- 2*upper
    lower <- 2*lower - 1
    add.rows <- (sol$model[,1]==i)*(sol$model[,2]>=lower)*(sol$model[,2]<=upper)
    rows <- c(rows, which(add.rows==1))
  }

  # remove the subtree subtended by the specified node and update the cluster assignment vector

  sol$Nodes <- sol$Nodes[-rows]

  sol$model <- sol$model[-rows,]

  sol$Parent <- sol$Parent[-rows]

  loci <- sol$model
  for(i in 1:max(sol$model[,1])){
    rows <- which(sol$model[,1]==i)
    loci[rows,2] <- rank(sol$model[rows,2])
  }

  for(i in 1:length(sol$Nodes)){
    sol$Nodes[[i]]$location <- loci[i,]
    sol$Nodes[[i]]$node <- sol$model[i,]
  }

  sol$cluster <- numeric(length(sol$cluster))

  for(i in 1:((length(sol$Nodes)+1)/2)) sol$cluster[sol$Nodes[[2*i-1]]$ixs] = i

  sol
}

### function tree_split allows a user to modify an existing clustering solution by further splitting
### a specified leaf node in an existing hierarchical model arising from one of the algorithms in the package.
## arguments:
# sol = clustering solution arising from one of the methods implemented in the package
# node = the node to be further split. Can be either the node number (based on the order of addition
#         of nodes to the model) or a vector describing the location of the node in the hierarchy
#         c(depth, position at depth)

## output is the updated clustering solution, which has the same form and type as the original solution

tree_split <- function(sol, node, ...){

  # store new parameter settings

  control <- list(...)

  # if location of node is given, find its node number

  if(length(node)>1){
    for(i in 1:length(sol$Nodes)){
      if(sum(sol$Nodes[[i]]$location==node)==2){
        node <- i
        break
      }
    }
    if(length(node)>1) stop('You must specify a node location within the given hierarchy')
  }

  # only leaf nodes can be further split. If an internal node is undesirable, and should be split
  # in a different way, then first prune the tree at that node, and then split it again using
  # split_leaf, with modified parameters

  is.leaf = 1-sum((sol$model[,1]==(sol$model[node,1]+1))*(sol$model[,2]==(2*sol$model[node,2])))
  if(!is.leaf) stop('only leaf nodes can be split within an existing hierarchy')

  # add the split at the specified node and find the details of the new child nodes. Then update
  # all of the structures contained in the solutio. Using the
  # same method as in the clustering algorithms themselves.

  if(sol$method=='MCDC'){

    # if new parameters were given, then modify the node being split before applying the partition

    if(length(control)>0){

      args_split <- sol$args

      if(!is.null(control$minsize)) args_split$minsize <- control$minsize
      if(!is.null(control$v0)) args_split$v0 <- control$v0
      if(!is.null(control$maxit)) args_split$maxit <- control$maxit
      if(!is.null(control$ftol)) args_split$ftol <- control$ftol

      c.split <- mch(sol$data[sol$Nodes[[node]]$ixs,], v0 = args_split$v0, minsize = args_split$minsize, maxit = args_split$maxit, ftol = args_split$ftol)

      sol$Nodes[[node]]$v <- c.split$v

      sol$Nodes[[node]]$b <- c.split$b

      sol$Nodes[[node]]$fval <- c.split$fval

      sol$Nodes[[node]]$params <- c.split$params

    }

    pass <- which(sol$data[sol$Nodes[[node]]$ixs,]%*%sol$Nodes[[node]]$v<sol$Nodes[[node]]$b)

    n.clust <- length(sol$Nodes)

    sol$Nodes[[n.clust+1]] <- list(ixs = sol$Nodes[[node]]$ixs[pass])

    sol$Nodes[[n.clust+2]] <- list(ixs = sol$Nodes[[node]]$ixs[-pass])

    sol$Parent <- c(sol$Parent, node, node)

    c.split <- mch(sol$data[sol$Nodes[[n.clust+1]]$ixs,], v0 = sol$args$v0, minsize = sol$args$minsize, maxit = sol$args$maxit, ftol = sol$args$ftol)

    sol$model <- rbind(sol$model, c(sol$model[node,1]+1, 2*sol$model[node,2]-1))

    sol$model <- rbind(sol$model, c(sol$model[node,1]+1, 2*sol$model[node,2]))

    sol$Nodes[[n.clust+1]]$v <- c.split$v

    sol$Nodes[[n.clust+1]]$b <- c.split$b

    sol$Nodes[[n.clust+1]]$fval <- c.split$fval

    sol$Nodes[[n.clust+1]]$params <- c.split$params

    sol$Nodes[[n.clust+1]]$node <- sol$model[n.clust+1,]

    c.split <- mch(sol$data[sol$Nodes[[n.clust+2]]$ixs,], v0 = sol$args$v0, minsize = sol$args$minsize, maxit = sol$args$maxit, ftol = sol$args$ftol)

    sol$Nodes[[n.clust+2]]$v <- c.split$v

    sol$Nodes[[n.clust+2]]$b <- c.split$b

    sol$Nodes[[n.clust+2]]$fval <- c.split$fval

    sol$Nodes[[n.clust+2]]$params <- c.split$params

    sol$Nodes[[n.clust+2]]$node <- sol$model[n.clust+2,]

    loci <- sol$model
    for(i in 1:max(sol$model[,1])){
      rows <- which(sol$model[,1]==i)
      loci[rows,2] <- rank(sol$model[rows,2])
    }

    sol$cluster <- numeric(length(sol$cluster))

    for(i in 1:((length(sol$Nodes)+1)/2)) sol$cluster[sol$Nodes[[2*i-1]]$ixs] = i

    for(i in 1:length(sol$Nodes)) sol$Nodes[[i]]$location <- loci[i,]

    return(sol)

  }
  else if(sol$method=='MDH'){

    # if new parameters were given, then modify the node being split before applying the partition

    if(length(control)>0){

      args_split <- sol$args

      if(!is.null(control$minsize)) args_split$minsize <- control$minsize
      if(!is.null(control$v0)) args_split$v0 <- control$v0
      if(!is.null(control$bandwidth)) args_split$bandwidth <- control$bandwidth
      if(!is.null(control$alphamin)) args_split$alphamin <- control$alphamin
      if(!is.null(control$alphamax)) args_split$alphamax <- control$alphamax
      if(!is.null(control$maxit)) args_split$maxit <- control$maxit
      if(!is.null(control$ftol)) args_split$ftol <- control$ftol

      c.split <- mdh(sol$data[sol$Nodes[[node]]$ixs,], v0 = args_split$v0, minsize = args_split$minsize, bandwidth = args_split$bandwidth, alphamin = args_split$alphamin, alphamax = args_split$alphamax, maxit = args_split$maxit, ftol = args_split$ftol)

      sol$Nodes[[node]]$v <- c.split$v

      sol$Nodes[[node]]$b <- c.split$b

      sol$Nodes[[node]]$fval <- c.split$fval

      sol$Nodes[[node]]$params <- c.split$params

      sol$Nodes[[node]]$rel.dep <- c.split$rel.dep

    }

    pass <- which(sol$data[sol$Nodes[[node]]$ixs,]%*%sol$Nodes[[node]]$v<sol$Nodes[[node]]$b)

    params = list()

    n <- nrow(sol$data)
    d <- ncol(sol$data)

    n.clust <- length(sol$Nodes)

    sol$Nodes[[n.clust+1]] <- list(ixs = sol$Nodes[[node]]$ixs[pass])

    sol$Nodes[[n.clust+2]] <- list(ixs = sol$Nodes[[node]]$ixs[-pass])

    sol$Parent <- c(sol$Parent, node, node)

    c.split <- mdh(sol$data[sol$Nodes[[n.clust+1]]$ixs,], v0 = sol$args$v0, minsize = sol$args$minsize, bandwidth = sol$args$bandwidth, alphamin = sol$args$alphamin, alphamax = sol$args$alphamax, maxit = sol$args$maxit, ftol = sol$args$ftol)

    sol$model <- rbind(sol$model, c(sol$model[node,1]+1, 2*sol$model[node,2]-1))

    sol$model <- rbind(sol$model, c(sol$model[node,1]+1, 2*sol$model[node,2]))

    sol$Nodes[[n.clust+1]]$v <- c.split$v

    sol$Nodes[[n.clust+1]]$b <- c.split$b

    sol$Nodes[[n.clust+1]]$fval <- c.split$fval

    sol$Nodes[[n.clust+1]]$params <- c.split$params

    sol$Nodes[[n.clust+1]]$rel.dep <- c.split$rel.dep

    sol$Nodes[[n.clust+1]]$node <- sol$model[n.clust+1,]

    c.split <- mdh(sol$data[sol$Nodes[[n.clust+2]]$ixs,], v0 = sol$args$v0, minsize = sol$args$minsize, bandwidth = sol$args$bandwidth, alphamin = sol$args$alphamin, alphamax = sol$args$alphamax, maxit = sol$args$maxit, ftol = sol$args$ftol)

    sol$Nodes[[n.clust+2]]$v <- c.split$v

    sol$Nodes[[n.clust+2]]$b <- c.split$b

    sol$Nodes[[n.clust+2]]$fval <- c.split$fval

    sol$Nodes[[n.clust+2]]$params <- c.split$params

    sol$Nodes[[n.clust+2]]$rel.dep <- c.split$rel.dep

    sol$Nodes[[n.clust+2]]$node <- sol$model[n.clust+2,]

    loci <- sol$model
    for(i in 1:max(sol$model[,1])){
      rows <- which(sol$model[,1]==i)
      loci[rows,2] <- rank(sol$model[rows,2])
    }

    sol$cluster <- numeric(length(sol$cluster))

    for(i in 1:((length(sol$Nodes)+1)/2)) sol$cluster[sol$Nodes[[2*i-1]]$ixs] = i

    for(i in 1:length(sol$Nodes)) sol$Nodes[[i]]$location <- loci[i,]

    return(sol)
  }
  else if(sol$method=='NCutH'){
    # if new parameters were given, then modify the node being split before applying the partition

    if(length(control)>0){

      args_split <- sol$args

      if(!is.null(control$minsize)) args_split$minsize <- control$minsize
      if(!is.null(control$v0)) args_split$v0 <- control$v0
      if(!is.null(control$s)) args_split$s <- control$s
      if(!is.null(control$maxit)) args_split$maxit <- control$maxit
      if(!is.null(control$ftol)) args_split$ftol <- control$ftol

      c.split <- ncuth(sol$data[sol$Nodes[[node]]$ixs,], v0 = args_split$v0, s = args_split$s, minsize = args_split$minsize, maxit = args_split$maxit, ftol = args_split$ftol)

      sol$Nodes[[node]]$v <- c.split$v

      sol$Nodes[[node]]$b <- c.split$b

      sol$Nodes[[node]]$fval <- c.split$fval

      sol$Nodes[[node]]$params <- c.split$params

    }


    pass <- which(sol$data[sol$Nodes[[node]]$ixs,]%*%sol$Nodes[[node]]$v<sol$Nodes[[node]]$b)

    n.clust <- length(sol$Nodes)

    sol$Nodes[[n.clust+1]] <- list(ixs = sol$Nodes[[node]]$ixs[pass])

    sol$Nodes[[n.clust+2]] <- list(ixs = sol$Nodes[[node]]$ixs[-pass])

    sol$Parent <- c(sol$Parent, node, node)

    c.split <- ncuth(sol$data[sol$Nodes[[n.clust+1]]$ixs,], v0 = sol$args$v0, s = sol$args$s, minsize = sol$args$minsize, maxit = sol$args$maxit, ftol = sol$args$ftol)

    sol$model <- rbind(sol$model, c(sol$model[node,1]+1, 2*sol$model[node,2]-1))

    sol$model <- rbind(sol$model, c(sol$model[node,1]+1, 2*sol$model[node,2]))

    sol$Nodes[[n.clust+1]]$v <- c.split$v

    sol$Nodes[[n.clust+1]]$b <- c.split$b

    sol$Nodes[[n.clust+1]]$fval <- c.split$fval

    sol$Nodes[[n.clust+1]]$params <- c.split$params

    sol$Nodes[[n.clust+1]]$node <- sol$model[n.clust+1,]

    c.split <- ncuth(sol$data[sol$Nodes[[n.clust+2]]$ixs,], v0 = sol$args$v0, s = sol$args$s, minsize = sol$args$minsize, maxit = sol$args$maxit, ftol = sol$args$ftol)

    sol$Nodes[[n.clust+2]]$v <- c.split$v

    sol$Nodes[[n.clust+2]]$b <- c.split$b

    sol$Nodes[[n.clust+2]]$fval <- c.split$fval

    sol$Nodes[[n.clust+2]]$params <- c.split$params

    sol$Nodes[[n.clust+2]]$node <- sol$model[n.clust+2,]

    loci <- sol$model
    for(i in 1:max(sol$model[,1])){
      rows <- which(sol$model[,1]==i)
      loci[rows,2] <- rank(sol$model[rows,2])
    }

    sol$cluster <- numeric(length(sol$cluster))

    for(i in 1:((length(sol$Nodes)+1)/2)) sol$cluster[sol$Nodes[[2*i-1]]$ixs] = i

    for(i in 1:length(sol$Nodes)) sol$Nodes[[i]]$location <- loci[i,]

    return(sol)
  }
}

