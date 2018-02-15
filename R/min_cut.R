#### Implements the NCutH method of D. Hofmeyr (2016), "Clustering by minimum cut hyperplanes," IEEE TPAMI.
#### The algorithm generates a divisive hierarchical clustering model using minimum normalised cut hyperplanes.
#### Differences in terminology from the above paper (for consistency with the package) are that here
#### each hyperplane is referred to as a NCutH, and the function to find such a hyperplane is ncuth.
#### The complete divisive clustering algorithm is herein referred to as NCutDC (function ncutdc), for
#### Normalised Cut Divisive Clustering.

### function ncutdc generates a complete clustering model using minimum normalised cut hyperplanes.
## arguments:
# X = dataset (matrix). each row is a datum. required
# K = number of clusters to extract (integer). required
# split.index = determines the order in which clusters are split (in decreasing
#       order of splitting indices). can be a function(v, X, P) of projection
#       vector v, data matrix X and list of parameters P. can also be one of
#       "size" (split the largest cluster), or "fval" (split the cluster with
#       the minimum normalised cut value). optional, default is "fval"
# v0 = initial projection direction(s). can be a matrix
#       in which each column is an initialisation to try.
#       can be a function of the data matrix (or subset
#       thereof corresponding to the cluster being split) which returns
#       a matrix in which each column is an initialisation.
#       optional, default is the first principal component
# s = used to compute scale parameter. a function f(X) of the data (or subset thereof) being split which returns a numeric in which
#      case scale = f(X). default is scale = 100*eigen(cov(X))$values[1]^.5*nrow(X)^(-0.2).
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
# maxit = maximum number of BFGS iterations for each value of alpha. optional, default is 50
# ftol = tolerance level for function value improvements in BFGS. optional, default is 1e-8

## output is a named list containing
# $cluster = cluster assignment vector
# $model = matrix containing the would-be location of each node (depth and position at depth) within a complete tree
# $Nodes = the clustering model. unnamed list each element of which is a named list containing the details of the associated node
# $data = the data matrix passed to mddc()
# $method = "NCutH" (used for plotting and model modification functions)
# $args = list of (functional) arguments passed to ncutdc

ncutdc <- function(X, K, split.index = NULL, v0 = NULL, s = NULL, minsize = NULL, verb = NULL, labels = NULL, maxit = NULL, ftol = NULL){
  
  if(is.data.frame(X)) X <- as.matrix(X)
  
  if(is.null(verb)) verb = 0

  # set parameters for clustering and optimisation

  if(is.null(split.index) || split.index=="fval") split.index <- function(v, X, P){
    if(nrow(rbind(c(), X))<=2) return(-Inf)
    1/(1+f_ncut(v, X, P))
  }
  else if(split.index=="size") split.index <- function(v, X, P) nrow(X)
  else if(!is.function(split.index)) stop('split.index must be one of "size" or "fval" or, a function of projection vector, data matrix and parameter list P with elements P$s and P$nmin')


  n <- nrow(X)
  d <- ncol(X)

  # obtain clusters and return cluster assignment vector

  split_indices <- numeric(2*K-1) - Inf

  # ixs represents a list of cluster indices (for each node in the hierarchy model structure)
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

  # stores the normalised cut values on each hyperplane

  ncuts <- numeric(2*K-1)

  # find the first hyperplane. if multiple initialisations were used, select the one with the minimum ncut value

  c.split <- ncuth(X, v0, s, minsize, verb, labels, maxit, ftol)
  ix.opt <- which.min(unlist(lapply(c.split, function(sol) sol$fval)))
  c.split <- c.split[[ix.opt]]

  # store the results in the above discussed objects

  split_indices[1] <- split.index(c.split$v, X, c.split$params)

  pass <- list(which(c.split$cluster==2))

  vs[1,] <- c.split$v

  bs[1] <- c.split$b

  pars[[1]] <- c.split$params

  ncuts[1] <- c.split$fval

  # repeat the divisive procedure until the correct number of clusters results

  while(length(ixs)<(2*K-1)){

    # split the leaf with the maximum split index

    id <- which.max(split_indices)
    split_indices[id] <- -Inf

    n.clust <- length(ixs)

    ixs[[n.clust+1]] <- ixs[[id]][pass[[id]]]

    ixs[[n.clust+2]] <- ixs[[id]][-pass[[id]]]

    c.split <- ncuth(X[ixs[[n.clust+1]],], v0, s, minsize, verb, labels[ixs[[n.clust+1]]], maxit, ftol)
    ix.opt <- which.min(unlist(lapply(c.split, function(sol) sol$fval)))
    c.split <- c.split[[ix.opt]]

    split_indices[n.clust+1] <- split.index(c.split$v, X[ixs[[n.clust+1]],], c.split$params)

    pass[[n.clust+1]] <- which(c.split$cluster==2)

    tree[n.clust+1,] <- c(tree[id,1] + 1, 2*tree[id,2]-1)

    vs[n.clust+1,] <- c.split$v

    bs[n.clust+1] <- c.split$b

    ncuts[n.clust+1] <- c.split$fval

    pars[[n.clust+1]] <- c.split$params
    
    Parent[n.clust+1] <- id

    c.split <- ncuth(X[ixs[[n.clust+2]],], v0, s, minsize, verb, labels[ixs[[n.clust+2]]], maxit, ftol)
    ix.opt <- which.min(unlist(lapply(c.split, function(sol) sol$fval)))
    c.split <- c.split[[ix.opt]]

    split_indices[n.clust+2] <- split.index(c.split$v, X[ixs[[n.clust+2]],], c.split$params)

    pass[[n.clust+2]] <- which(c.split$cluster==2)

    tree[n.clust+2,] <- c(tree[id,1] + 1, 2*tree[id,2])

    vs[n.clust+2,] <- c.split$v

    bs[n.clust+2] <- c.split$b

    ncuts[n.clust+2] <- c.split$fval

    pars[[n.clust+2]] <- c.split$params
    
    Parent[n.clust+2] <- id
  }

  # create cluster assignment vector for the complete hierarchy

  asgn <- numeric(n)+1
  for(i in 1:(K-1)) asgn[ixs[[2*i]]] <- i+1

  # determine actual locations of the nodes (not their would-be location within a complete tree)

  loci <- tree
  for(i in 1:max(tree[,1])){
    rows <- which(tree[,1]==i)
    loci[rows,2] <- rank(tree[rows,2])
  }

  # create the detailed cluster hierarchy

  Nodes <- list()
  for(i in 1:length(ixs)) Nodes[[i]] <- list(ixs = ixs[[i]], v = vs[i,], b = bs[i], params = pars[[i]], fval = ncuts[i], node = tree[i,], location = loci[i,])

  list(cluster = asgn, model = tree, Parent = Parent, Nodes = Nodes, data = X, method = 'NCutH', args = list(v0 = v0, s = s, split.index = split.index, minsize = minsize, maxit = maxit, ftol = ftol))
}



### function ncuth() finds minimum normalised cut hyperplanes
## arguments:
# X = dataset (matrix). each row is a datum. required
# v0 = initial projection direction(s). can be a matrix
#       in which each column is an initialisation to try.
#       can be a function of the data matrix (or subset
#       thereof corresponding to the cluster being split) which returns
#       a matrix in which each column is an initialisation.
#       optional, default is the first principal component
# s = positive numeric scaling parameter (sigma). optional, default is s = 100*eigen(cov(X))$values[1]^.5*nrow(X)^(-0.2)
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
# fval = the normalised cut across H(v, b)
# params = list of parameters used to find H(v, b)


ncuth <- function(X, v0 = NULL, s = NULL, minsize = NULL, verb = NULL, labels = NULL, maxit = NULL, ftol = NULL){

  if(is.data.frame(X)) X <- as.matrix(X)
  
  params = list()

  if(is.null(minsize)) params$nmin <- 1
  else if(is.function(minsize)) params$nmin <- minsize(X)
  else if(is.numeric(minsize) && length(minsize)==1) params$nmin <- minsize
  else stop('minsize must be a positive integer or a function of the data being split')

  if(is.vector(X)) n <- 1
  else n <- nrow(X)

  # if the data contain fewer than 2*nmin points then don't split

  if(n<(2*params$nmin)) return(list(list(cluster = numeric(n)+1, v = numeric(ncol(rbind(c(),X)))+1, b = 0, params = list(s = 100, nmin = params$nmin), fval = Inf, method = 'NCutH')))

  
  # if labels are supplied, ensure they are integers 1:K (K the number of classes)
  
  if(!is.null(labels)){
    lab_new <- numeric(length(labels))
    u <- unique(labels)
    for(i in 1:length(u)) lab_new[which(labels==u[i])] = i
    labels <- lab_new
  }
  
  # set up parameters for optimisation

  if(is.null(verb)) verb = 0

  if(is.null(maxit)) maxit <- 50

  if(is.null(ftol)) ftol <- 1e-8

  if(is.null(s)){
    if(ncol(X)>2) params$s <- rARPACK::eigs_sym(cov(X), 1)$values[1]^.5*100/nrow(X)^.2
    else params$s <- eigen(cov(X))$values[1]^.5*100/nrow(X)^.2
  }
  else if(is.numeric(s)) params$s <- s
  else if(is.function(s)) params$s <- s(X)
  else stop('s must be numeric or a function of the data being split')

  if(is.null(v0)){
    if(ncol(X)>2) E <- matrix(rARPACK::eigs_sym(cov(X), 1)$vectors, ncol = 1)
    else E <- matrix(eigen(cov(X))$vectors[,1], ncol = 1)
  }
  else if (is.function(v0)) E <- cbind(c(), v0(X))
  else if (is.vector(v0)) E <- matrix(v0, ncol = 1)
  else E <- cbind(c(), v0)

  # store details of the hyperplanes arising from each initialisation

  output <- list()

  for(i in 1:ncol(E)){
    v <- ncutpp(E[,i], X, params, verb, labels, maxit, ftol)

    b <- ncut_b(v, X, params)

    fval <- f_ncut(v, X, params)

    pass <- X%*%v<b

    output[[i]] <- list(cluster = pass+1, v = v, b = b, params = params, fval = fval, method = 'NCutH')

  }

  output

}

### function f_ncut evaluates the projection index for ncuth
## arguments:
# v = projection vector
# X = data matrix
# P = list of parameters containing (at least)
# $s = scaling constant
# $nmin = minimum cluster size

## output is a scalar, the value of the normalised cut across the best hyperplane orthogonal to v

f_ncut <- function(v, X, P){

  # project the data onto v and sort in increasing order

  x <- X%*%v/norm_vec(v)/P$s
  srt <- sort(x)
  n <- length(x)

  # determine the normalised cut at each projected point and return the minimum (see TPAMI paper for details)

  EG <- exp(srt[1]-srt)
  CSEG <- cumsum(EG)
  CSiEG <- cumsum(1/EG)
  Cut <- (CSEG[n]-CSEG[P$nmin:(n-P$nmin)])*CSiEG[P$nmin:(n-P$nmin)]
  Deg <- ((cumsum(EG*CSiEG) + cumsum((CSEG[n]-CSEG)/EG)))
  min((Cut*(1/Deg[P$nmin:(n-P$nmin)]+1/(Deg[n]-Deg[P$nmin:(n-P$nmin)]))))
}

### function df_ncut evaluates the gradient of the projection index for ncuth, provided the
### optimal hyperplane orthogonal to v is unique and the projected point at the optimum is unique
## arguments:
# v = projection vector
# X = data matrix
# P = list of parameters containing (at least)
# $s = scaling constant
# $nmin = minimum cluster size

## output is a vector. the gradient of the projection index at v

df_ncut <- function(v, X, P){

  # project the data onto v and determine their ordering

  x <- X%*%v/norm_vec(v)/P$s
  o <- order(x)
  srt <- x[o]
  n <- length(x)

  # find the location of the optimal hyperplane orthogonal to v

  EG <- exp(srt[1]-srt)
  CSEG <- cumsum(EG)
  CSiEG <- cumsum(1/EG)
  Cut <- (CSEG[n]-CSEG[P$nmin:(n-P$nmin)])*CSiEG[P$nmin:(n-P$nmin)]
  Deg <- ((cumsum(EG*CSiEG) + cumsum((CSEG[n]-CSEG)/EG)))
  w <- which.min((Cut*(1/Deg[P$nmin:(n-P$nmin)]+1/(Deg[n]-Deg[P$nmin:(n-P$nmin)])))) + P$nmin - 1

  # compute the gradient of the normalised cut (see TPAMI paper for details)

  Ek <- exp(srt-srt[w])
  Eki <- 1/Ek
  CSEk <- cumsum(Ek)
  CSEki <- cumsum(Eki)
  ds <- numeric(n)
  if(w>2){
    ds[1] <- (CSEki[n]-CSEki[w])*Ek[1]*(1/Deg[w]+1/(Deg[n]-Deg[w]))-Cut[w-P$nmin+1]*((2*Ek[1]*(CSEki[w]-CSEki[1])+Ek[1]*(CSEki[n]-CSEki[w]))/Deg[w]^2+Ek[1]*(CSEki[n]-CSEki[w])/(Deg[n]-Deg[w])^2)
    ds[2:(w-1)] <- (CSEki[n]-CSEki[w])*Ek[2:(w-1)]*(1/Deg[w]+1/(Deg[n]-Deg[w]))-Cut[w-P$nmin+1]*((2*Ek[2:(w-1)]*(CSEki[w]-CSEki[2:(w-1)])-2*Eki[2:(w-1)]*CSEk[1:(w-2)]+Ek[2:(w-1)]*(CSEki[n]-CSEki[w]))/Deg[w]^2+Ek[2:(w-1)]*(CSEki[n]-CSEki[w])/(Deg[n]-Deg[w])^2)
    ds[w] <- (CSEki[n]-CSEki[w])*(1/Deg[w]+1/(Deg[n]-Deg[w]))-Cut[w-P$nmin+1]*((CSEki[n]-CSEki[w]-2*CSEk[w-1])/Deg[w]^2+(CSEki[n]-CSEki[w])/(Deg[n]-Deg[w])^2)
    ds[(w+1):n] <- -CSEk[w]*Eki[(w+1):n]*(1/Deg[w]+1/(Deg[n]-Deg[w]))-Cut[w-P$nmin+1]*(-Eki[(w+1):n]*CSEk[w]/Deg[w]^2+(2*Ek[(w+1):n]*(CSEki[n]-CSEki[(w+1):n])-2*Eki[(w+1):n]*(CSEk[w:(n-1)]-CSEk[w])-Eki[(w+1):n]*CSEk[w])/(Deg[n]-Deg[w])^2)
  }
  else{
    ds <- sapply(1:n, function(l){
      if(l==1){
        d1 <- (CSEki[n]-CSEki[w])*Ek[l]
        d2 <- 2*Ek[l]*(CSEki[w]-CSEki[l])+Ek[l]*(CSEki[n]-CSEki[w])
        d3 <- Ek[l]*(CSEki[n]-CSEki[w])
        d1*(1/Deg[w]+1/(Deg[n]-Deg[w]))-Cut[w-P$nmin+1]*(d2/Deg[w]^2+d3/(Deg[n]-Deg[w])^2)
      }
      else if(l<w){
        d1 <- (CSEki[n]-CSEki[w])*Ek[l]
        d2 <- 2*Ek[l]*(CSEki[w]-CSEki[l])-2*Eki[l]*CSEk[l-1]+Ek[l]*(CSEki[n]-CSEki[w])
        d3 <- Ek[l]*(CSEki[n]-CSEki[w])
        d1*(1/Deg[w]+1/(Deg[n]-Deg[w]))-Cut[w-P$nmin+1]*(d2/Deg[w]^2+d3/(Deg[n]-Deg[w])^2)
      }
      else if(l>w){
        d1 <- -CSEk[w]*Eki[l]
        d2 <- -Eki[l]*CSEk[w]
        d3 <- 2*Ek[l]*(CSEki[n]-CSEki[l])-2*Eki[l]*(CSEk[l-1]-CSEk[w])-Eki[l]*CSEk[w]
        d1*(1/Deg[w]+1/(Deg[n]-Deg[w]))-Cut[w-P$nmin+1]*(d2/Deg[w]^2+d3/(Deg[n]-Deg[w])^2)
      }
      else{
        d1 <- CSEki[n]-CSEki[w]
        d2 <- CSEki[n]-CSEki[w]-2*CSEk[w-1]
        d3 <- CSEki[n] - CSEki[w]
        d1*(1/Deg[w]+1/(Deg[n]-Deg[w]))-Cut[w-P$nmin+1]*(d2/Deg[w]^2+d3/(Deg[n]-Deg[w])^2)
      }
    })
  }
  nv <- norm_vec(v)
  dv <- (X[o,]/nv-((X[o,])%*%v)%*%t(v)/nv^3)/P$s
  ds%*%dv
}

### function ncut_b determines the location of the optimal hyperplane orthogonal to projection vector v
### That is, the value of b making H(v, b) an optimal hyperplane
## arguments
# v = projection vector
# X = data matrix
# P = list of parameters containing (at least)
# $s = scaling constant
# $nmin = minimum cluster size

## output is a scalar. the location of the optimal hyperplane orthogonal to v

ncut_b <- function(v, X, P){

  # functionality follows closely that of f_ncut

  x <- X%*%v/norm_vec(v)/P$s
  srt <- sort(x)
  n <- length(x)
  EG <- exp(srt[1]-srt)
  CSEG <- cumsum(EG)
  CSiEG <- cumsum(1/EG)
  Cut <- (CSEG[n]-CSEG[P$nmin:(n-P$nmin)])*CSiEG[P$nmin:(n-P$nmin)]
  Deg <- ((cumsum(EG*CSiEG) + cumsum((CSEG[n]-CSEG)/EG)))
  w <- which.min((Cut*(1/Deg[P$nmin:(n-P$nmin)]+1/(Deg[n]-Deg[P$nmin:(n-P$nmin)]))))+P$nmin-1
  (srt[w] + srt[w+1])/2*P$s
}

### function ncutpp performs projection pursuit based on the normalised cut objective to find the optimal
### solution. Acts as a gateway to ppclust.optim, providing appropriate arguments for the ncut objective
## arguments:
# v = initial projection vector
# X = data matrix
# P = list of parameters containing (at least)
# $s = scaling parameter
# $nmin = minimum cluster size
# verb = verbosity level. see details in paper or at function ncutdc or ncuth
# labels = vector of class labels used only in plotting for verb> 0
# maxit = maximum number of iterations for optimisation
# ftol = relative tolerance level for the location of an optimum

ncutpp <- function(v, X, P, verb, labels, maxit, ftol){
  v_new <- ppclust.optim(v, f_ncut, df_ncut, X, P, ncut_b, verbosity = verb, labels = labels, method = 'NCutH', maxit = maxit, ftol = ftol)$par
  v_new/norm_vec(v_new)
}


