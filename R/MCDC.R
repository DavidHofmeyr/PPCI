#### Implements the Maximum Clusterability Divisive Clustering (MCDC) algorithm of
#### Hofmeyr and Pavlidis (2015), IEEE SSCI CIDM

### function mcdc generates a divisive hierarchical clustering model using hyperplanes which maximise the
### variance ratio clusterability measure across them
## arguments:
# X = dataset (matrix). each row is a datum. required
# K = number of clusters to extract (integer). required
# split.index = determines the order in which clusters are split (in decreasing
#       order of splitting indices). can be a function(v, X, P) of projection
#       vector v, data matrix X and list of parameters P. can also be one of
#       "size" (split the largest cluster), "fval" (split the cluster with
#       the maximum variance ratio value), or "Fdist" (indices determined by the non-central
#       F-distribution. See SSCI paper for details. slight difference from the paper is that
#       when the data size is above 2000 cluster size is used instead. This is because the naive
#       estimation of the model degrees of freedom has been seen to be unreliable when the number
#       of data is large). optional, default is "Fdist"
# v0 = initial projection direction(s). can be a matrix
#       in which each column is an initialisation to try.
#       can be a function of the data matrix (or subset
#       thereof corresponding to the cluster being split) which returns
#       a matrix in which each column is an initialisation.
#       optional, default is the vector joining the means of a 2-means solution
# minsize = minimum cluster size. Can be either integer or a function f(X) returning an integer. Default is 1. Throughout the projection pursuit no cuts which result in a cluster smaller
#            than minsize are allowed. This is achieved by considering only partitions in (v%*%X)[minsize:(n-minsize+1)].
# verb = verbosity level. verb == 0 produces no output. verb == 1 produces plots of the
#         projected data during each optimisation. verb == 2 adds to these plots information
#         about the function value, and quality of split (if labels are supplied).
#         verb == 3 creates a folder in the working directory and saves all plots produced for verb == 2.
#         optional, default is 3
# labels = vector of class labels. Only used for producing plots, not in the allocation of
#         data to clusters. optional, default is NULL (plots do not indicate true class membership
#         when true labels are unknown)
# maxit = maximum number of BFGS iterations for each value of alpha. optional, default is 15
# ftol = tolerance level for function value improvements in BFGS. optional, default is 1e-5

## output is a named list containing
# $cluster = cluster assignment vector
# $model = matrix containing the would-be location of each node (depth and position at depth) within a complete tree
# $Nodes = the clustering model. unnamed list each element of which is a named list containing the details of the associated node
# $data = the data matrix passed to mcdc()
# $method = "MCDC" (used for plotting and model modification functions)
# $args = list of (functional) arguments passed to ncutdc

mcdc <- function(X, K, v0 = NULL, split.index = NULL, minsize = NULL, verb = NULL, labels = NULL, maxit = NULL, ftol = NULL){

  if(is.data.frame(X)) X <- as.matrix(X)
  
  if(is.null(verb)) verb = 0

  # set parameters for clustering and optimisation

  n <- nrow(X)
  d <- ncol(X)

  if(is.null(split.index) || split.index=='Fdist'){
    if(n<2000){
      split.index <- function(v, X, P){
        if(nrow(rbind(c(), X))<=2) return(-Inf)
        VR <- f_mc(v, X, P)
        n <-nrow(X)
        d <- ncol(X)
        al <- min(n, d+1)
        beta <- max(0, n-d-1)
        if(beta==0) return(0)
        pf(VR*beta/al, al, beta, ncp = n) + 1e-30*sqrt(n)*VR
      }
    }
    else split.index <- function(v, X, P) nrow(X)
  }
  else if(split.index=='size') split.index <- function(v, X, P) nrow(X)
  else if(split.index=='fval') split.index <- function(v, X, P) f_mc(v, X, P)
  else if(!is.function(split.index)) stop('split.index must be a function of projection vector, data matrix and parameter list P with elements P$nmin')

  # obtain clusters and cluster hierarchy

  # split_indices used to select the order to partition nodes/clusters

  split_indices <- numeric(2*K-1) - Inf

  # ixs stores the data associated with each node in the model

  ixs <- list(1:n)

  # tree stores the location (depth, breadth) in the model of each node

  tree <- matrix(0, (2*K-1), 2)
  tree[1,] <- c(1, 1)

  # Parent stores the parent node number of each node (The parent of the root node is 0)

  Parent <- numeric(2*K-1)

  # stores hyperplane separators for each node, v and b

  vs <- matrix(0, (2*K-1), d)
  bs <- numeric(2*K-1)

  # stores the parameters used in each optimisation

  pars <- list()
  VRS <- numeric(2*K-1)

  # determine the optimal hyperplane(s) at the root node and select that with the maximum variance ratio

  c.split <- mch(X, v0, minsize, verb, labels, maxit, ftol)
  ix.opt <- which.max(unlist(lapply(c.split, function(sol) sol$fval)))
  c.split <- c.split[[ix.opt]]

  # store the results in the above discussed objects

  split_indices[1] <- split.index(c.split$v, X, c.split$params)

  pass <- list(which(c.split$cluster==2))

  vs[1,] <- c.split$v

  bs[1] <- c.split$b

  pars[[1]] <- c.split$params

  VRS[1] <- c.split$fval

  # repeatedly apply binary partitions until the desired number of clusters results

  while(length(ixs)<(2*K-1)){

    # select the leaf with the greatest split index

    id <- which.max(split_indices)
    split_indices[id] <- -Inf

    n.clust <- length(ixs)

    ixs[[n.clust+1]] <- ixs[[id]][pass[[id]]]

    ixs[[n.clust+2]] <- ixs[[id]][-pass[[id]]]

    c.split <- mch(X[ixs[[n.clust+1]],], v0, minsize, verb, labels[ixs[[n.clust+1]]], maxit, ftol)
    ix.opt <- which.max(unlist(lapply(c.split, function(sol) sol$fval)))
    c.split <- c.split[[ix.opt]]

    split_indices[n.clust+1] <- split.index(c.split$v, X[ixs[[n.clust+1]],], c.split$params)

    pass[[n.clust+1]] <- which(c.split$cluster==2)

    tree[n.clust+1,] <- c(tree[id,1] + 1, 2*tree[id,2]-1)

    vs[n.clust+1,] <- c.split$v

    bs[n.clust+1] <- c.split$b

    VRS[n.clust+1] <- c.split$fval

    pars[[n.clust+1]] <- c.split$params

    Parent[n.clust+1] <- id

    c.split <- mch(X[ixs[[n.clust+2]],], v0, minsize, verb, labels[ixs[[n.clust+2]]], maxit, ftol)
    ix.opt <- which.max(unlist(lapply(c.split, function(sol) sol$fval)))
    c.split <- c.split[[ix.opt]]

    split_indices[n.clust+2] <- split.index(c.split$v, X[ixs[[n.clust+2]],], c.split$params)

    pass[[n.clust+2]] <- which(c.split$cluster==2)

    tree[n.clust+2,] <- c(tree[id,1] + 1, 2*tree[id,2])

    vs[n.clust+2,] <- c.split$v

    bs[n.clust+2] <- c.split$b

    VRS[n.clust+2] <- c.split$fval

    pars[[n.clust+2]] <- c.split$params

    Parent[n.clust+2] <- id
  }

  # determine cluster assignment vector

  asgn <- numeric(n) + 1
  for(i in 1:(K-1)) asgn[ixs[[2*i]]] <- i+1

  # find the actual location of each node in the hierarchy

  loci <- tree
  for(i in 1:max(tree[,1])){
    rows <- which(tree[,1]==i)
    loci[rows,2] <- rank(tree[rows,2])
  }

  # store the details of all hyperplanes used in the hierarchical model

  Nodes <- list()
  for(i in 1:length(ixs)) Nodes[[i]] <- list(ixs = ixs[[i]], v = vs[i,], b = bs[i], params = pars[[i]], fval = VRS[i], node = tree[i,], location = loci[i,])

  list(cluster = asgn, model = tree, Parent = Parent, Nodes = Nodes, data = X, method = 'MCDC', args = list(v0 = v0, split.index = split.index, minsize = minsize, maxit = maxit, ftol = ftol))
}

### function f_mc evaluates the projection index for mcdc
## arguments:
# v = projection vector
# X = data matrix
# P = list of parameters containing (at least)
# $nmin = minimum cluster size

## output is a scalar, the variance ratio clusterability of the optimal partition by a hyperplane orthogonal to v

f_mc <- function(v, X, P){

  # compute the projected points and sort them in increasing order

  p <- sort(X%*%v/norm_vec(v))
  CS <- cumsum(p)
  n <- length(p)

  # find the variance ratio at each point and return the maximum

  V <- sum((p-CS[n]/n)^2)/(n-1)
  ixs <- P$nmin:(n-P$nmin)
  bc <- (n-ixs)/ixs*((CS[n]-CS[ixs])/(n-ixs)-CS[n]/n)^2
  max(bc/(V-bc))
}

### function df_mc evaluates the gradient of the projection index for mcdc
### the gradient is valid when the optimal partition is unique and the projected
### point at the optimum is unique. This is a.e. w.r.t the lebesuge measure
## arguments:
# v = projection vector
# X = data matrix
# P = list of parameters containing (at least)
# $nmin = minimum cluster size

## output is a vector, the gradient of the variance ratio of the optimal hyperplane
## orthogonal to v

df_mc <- function(v, X, P){

  # compute the projected points and their ordering

  p <- X%*%v/norm_vec(v)
  o <- order(p)
  CS <- cumsum(p[o])
  n <- length(p)

  # determine the location of the optimal split

  V <- sum((p-CS[n]/n)^2)/(n-1)
  ixs <- P$nmin:(n-P$nmin)
  bc <- (n-ixs)/ixs*((CS[n]-CS[ixs])/(n-ixs)-CS[n]/n)^2
  VRS <- bc/(V-bc)
  b <- p[o][ixs][which.max(VRS)]

  # compute the gradient (which differs depending which side of the split each projected point lies, ix1 vs ix2)

  ix1 <- which(p<=b)
  ix2 <- which(p>b)
  m1 <- mean(p[ix1])
  m2 <- mean(p[ix2])
  m <- mean(p)
  bc <- length(ix1)/n*(m1-m)^2+length(ix2)/n*(m2-m)^2
  dp1 <- (bc*(2/n*(m-m1)-2/n*(p[ix1]-m1))+2/n*(m1-m)*V)/(V-bc)^2
  dp2 <- (bc*(2/n*(m-m2)-2/n*(p[ix2]-m2))+2/n*(m2-m)*V)/(V-bc)^2
  dp <- numeric(n)
  dp[ix1] <- dp1
  dp[ix2] <- dp2
  nv <- norm_vec(v)
  dv <- (X/nv-((X)%*%v)%*%t(v)/nv^3)
  dp%*%dv
}


### function mc_b finds the location of the optimal hyperplane orthogonal to v. That is,
### the value of b which makes H(v, b) an optimal hyperplane

mc_b <- function(v, X, P){

  # follows essentially the same procedure as evaluating the projection index

  p <- X%*%v/norm_vec(v)
  o <- order(p)
  CS <- cumsum(p[o])
  n <- length(p)
  V <- sum((p-CS[n]/n)^2)/(n-1)
  ixs <- P$nmin:(n-P$nmin)
  bc <- (n-ixs)/ixs*((CS[n]-CS[ixs])/(n-ixs)-CS[n]/n)^2
  VRS <- bc/(V-bc)
  w <- which.max(VRS)
  (p[o][ixs][w] + p[o][ixs][w+1])/2
}

### function mcpp performs projection pursuit based on variance ratio objective. The function
### acts as a gateway to the optimisation function ppclust.optim, providing appropriate arguments
### for mcdc
## arguments:
# v = initial projection vector
# X = data matrix
# P = list of paramateters containing (at least)
# $nmin = minimum cluster size
# verb = verbosity level. See details in paper or at function mcdc/mch
# maxit = maximum number of iterations in optimisation
# ftol = relative tolerance level for convergence of gradient based optimisation

## output is the optimal projection vector

mcpp <- function(v, X, P, verb, labels, maxit, ftol){
  v <- ppclust.optim(v, f_mc, df_mc, X, P, mc_b, verbosity = verb, labels = labels, method = 'MCDC', maxit = maxit, ftol = ftol)$par
  return(v/norm_vec(v))
}


### function mch() finds maximum variance ratio hyperplanes
## arguments:
# X = dataset (matrix). each row is a datum. required
# v0 = initial projection direction(s). can be a matrix
#       in which each column is an initialisation to try.
#       can be a function of the data matrix (or subset
#       thereof corresponding to the cluster being split) which returns
#       a matrix in which each column is an initialisation.
#       optional, default is the vector joining the means of a 2-means clustering
# verb = verbosity level. verb == 0 produces no output. verb == 1 produces plots of the
#         projected data during each optimisation. verb == 2 adds to these plots information
#         about the function value, relative depth and quality of split (if labels are supplied).
#         verb == 3 creates a folder in the working directory and saves all plots produced for verb == 2.
#         optional, default is 3
# labels = vector of class labels. Only used for producing plots, not in the allocation of
#         data to clusters. optional, default is NULL (plots do not indicate true class membership
#         when true labels are unknown)
# maxit = maximum number of BFGS iterations for each value of alpha. optional, default is 15
# ftol = tolerance level for function value improvements in BFGS. optional, default is 1e-5

## output is a list of lists, the i-th stores the details of the optimal hyperplane
## arising from the initialisation at v0[,i]. Each element has contains
# $cluster = the cluster assignment vector
# $v = the optimal projection vector
# $b = the value of b making H(v, b) the optimal hyperplane
# fval = the variance ratio across H(v, b)
# params = list of parameters used to find H(v, b)

mch <- function(X, v0 = NULL, minsize = NULL, verb = NULL, labels = NULL, maxit = NULL, ftol = NULL){
  
  if(is.data.frame(X)) X <- as.matrix(X)
  
  params = list()

  if(is.null(minsize)) params$nmin <- 1
  else if(is.function(minsize)) params$nmin <- minsize(X)
  else if(is.numeric(minsize) && length(minsize)==1) params$nmin <- minsize
  else stop('minsize must be a positive integer or a function of the data being split')

  # if labels are supplied, ensure they are integers 1:K (K the number of classes)

  if(!is.null(labels)){
    lab_new <- numeric(length(labels))
    u <- unique(labels)
    for(i in 1:length(u)) lab_new[which(labels==u[i])] = i
    labels <- lab_new
  }

  # if there are fewer than 2*minsize data, do not split

  if(is.vector(X)) n <- 1
  else n <- nrow(X)
  if(n<(2*params$nmin)) return(list(list(cluster = c(0, numeric(n-1) + 1), v = numeric(ncol(rbind(c(), X))) + 1, b = 0, params = list(nmin=1), fval = -Inf)))

  if(is.null(verb)) verb = 0


  # set up parameters for optimisation

  if(is.null(maxit)) maxit <- 50

  if(is.null(ftol)) ftol <- 1e-8

  if(is.null(v0)){
    km <- kmeans(X, 2, nstart = 10)
    E <- cbind(c(), km$centers[1,]-km$centers[2,])
  }
  else if (is.function(v0)) E <- cbind(c(), v0(X))
  else E <- cbind(c(), v0)

  output <- list()

  for(i in 1:ncol(E)){
    v <- mcpp(E[,i], X, params, verb, labels, maxit, ftol)

    b <- mc_b(v, X, params)

    fval <- f_mc(v, X, params)

    pass <- X%*%v<b

    output[[i]] <- list(cluster = pass + 1, v = v, b = b, params = params, fval = fval, method = 'MCDC')
  }

  output

}

### function norm_vec computes the euclidean norm of a vector. This function is used by all methods in the package
## arguments:
# v = numeric vector

## output is a scalar, the euclidean norm of the vector

norm_vec <- function(v) sqrt(sum(v^2))
