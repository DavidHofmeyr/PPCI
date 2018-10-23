#### Implements the MDPP algorithm of Pavlidis et. al (2016) to find minimum density hyperplanes (MDHs).
#### In addition, multiple MDHs can be combined to produce a complete clustering using a divisive hierarchical
#### algorithm.



### function mddc() generates a complete clustering model using MDHs
## arguments:
# X = dataset (matrix). each row is a datum. required
# K = number of clusters to extract (integer). required
# split.index = determines the order in which clusters are split (in decreasing
#       order of splitting indices). can be a function(v, X, P) of projection
#       vector v, data matrix X and list of parameters P. can also be one of
#       "size" (split the largest cluster), "fval" (split the cluster with
#       the minimum density hyperplane), or "rdepth" (split the cluster with
#       the greatest relative depth). optional, default is "size"
# v0 = initial projection direction(s). can be a matrix
#       in which each column is an initialisation to try.
#       can be a function of the data matrix (or subset
#       thereof corresponding to the cluster being split) which returns
#       a matrix in which each column is an initialisation.
#       optional, default is the first principal component
# bandwidth = used to compute bandwidth parameter (h). a function f(X) of
#       the data (or subset) being pslit. optional, default is bandwidth(X) = 0.9*eigen(cov(X))$values[1]^.5/nrow(X)^.2
# alphamax = maximum width of constraint F(v) = [mu - alphamax*sd, mu+alphamax*sd]
#       optional, default is 0.9
# verb = verbosity level. verb == 0 produces no output. verb == 1 produces plots of the
#         projected data during each optimisation. verb == 2 adds to these plots information
#         about the function value, relative depth and quality of split (if labels are supplied).
#         verb == 3 creates a folder in the working directory and saves all plots produced for verb == 2.
#         optional, default is 0
# labels = vector of class labels. Only used for producing plots, not in the allocation of
#         data to clusters. optional, default is NULL (plots do not indicate true class membership
#         when true labels are unknown)
# maxit = maximum number of BFGS iterations for each value of alpha. optional, default is 15
# ftol = tolerance level for function value improvements in BFGS. optional, default is 1e-5

## output is a named list containing
# $cluster = cluster assignment vector
# $model = matrix containing the would-be location of each node (depth and position at depth) within a complete tree
# $Nodes = the clustering model. unnamed list each element of which is a named list containing the details of the associated node
# $data = the data matrix passed to mddc()
# $method = "MDH" (used for plotting and model modification functions)
# $args = list of (functional) arguments passed to mddc

mddc <- function(X, K, minsize = NULL, split.index = NULL, v0 = NULL, bandwidth = NULL, alphamin = NULL, alphamax = NULL, verb = NULL, labels = NULL, maxit = NULL, ftol = NULL){

  if(is.data.frame(X)) X <- as.matrix(X)

  if(is.null(verb)) verb = 0
  # set parameters for clustering and optimisation

  if(is.null(split.index)) split.index <- function(v, X, P) nrow(X)
  else if(is.character(split.index)){
    if(split.index=='size') split.index <- function(v, X, P) nrow(X)
    else if(split.index=='fval') split.index <- function(v, X, P) 1/f_md(v, X, P)
    else if(split.index=='rdepth') split.index <- function(v, X, P) md_reldepth(v, X, P)
    else stop('split.index must be one of "size", "fval", or "rdepth"')
  }
  else if(!is.function(split.index)) stop('split.index must be one of "size", "fval", or "rdepth" or a function f(v, X, P) of projection vector v, data X and list of parameters P')

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

  # store the integrated density on each hyperplane and the corresponding relative depth respectively
  dens <- numeric(2*K-1)
  rel.deps <- numeric(2*K-1)

  # compute the split(s) at the root node, and choose the solution with the maximum relative depth

  c.split <- mdh(X, v0, minsize, bandwidth, alphamin, alphamax, verb, labels, maxit, ftol)

  if(c.split$rel.dep>0) split_indices[1] <- split.index(c.split$v, X, c.split$params)
  else split_indices[1] <- -Inf

  pass <- list(which(c.split$cluster==2))

  vs[1,] <- c.split$v

  bs[1] <- c.split$b

  rel.deps[1] <- c.split$rel.dep

  dens[1] <- c.split$fval

  pars[[1]] <- c.split$params


  # repeat binary partitions until the desired number of clusters have been generated

  while(length(ixs)<(2*K-1) && max(split_indices) > -Inf){
    id <- which.max(split_indices)

    split_indices[id] <- -Inf

    n.clust <- length(ixs)

    ixs[[n.clust+1]] <- ixs[[id]][pass[[id]]]

    ixs[[n.clust+2]] <- ixs[[id]][-pass[[id]]]

    c.split <- mdh(X[ixs[[n.clust+1]],], v0, minsize, bandwidth, alphamin, alphamax, verb, labels[ixs[[n.clust+1]]], maxit, ftol)

    if(c.split$rel.dep==0) split_indices[n.clust+1] <- -Inf
    else split_indices[n.clust+1] <- split.index(c.split$v, X[ixs[[n.clust+1]],], c.split$params)

    pass[[n.clust+1]] <- which(c.split$cluster==2)

    vs[n.clust+1,] <- c.split$v

    bs[n.clust+1] <- c.split$b

    rel.deps[n.clust+1] <- c.split$rel.dep

    dens[n.clust+1] <- c.split$fval

    pars[[n.clust+1]] <- c.split$params

    tree[n.clust+1,] <- c(tree[id,1] + 1, 2*tree[id,2]-1)

    Parent[n.clust+1] <- id

    c.split <- mdh(X[ixs[[n.clust+2]],], v0, minsize, bandwidth, alphamin, alphamax, verb, labels[ixs[[n.clust+2]]], maxit, ftol)

    if(c.split$rel.dep==0) split_indices[n.clust+2] <- -Inf
    else split_indices[n.clust+2] <- split.index(c.split$v, X[ixs[[n.clust+2]],], c.split$params)

    pass[[n.clust+2]] <- which(c.split$cluster==2)

    vs[n.clust+2,] <-c.split$v

    bs[n.clust+2] <- c.split$b

    rel.deps[n.clust+2] <- c.split$rel.dep

    dens[n.clust+2] <- c.split$fval

    pars[[n.clust+2]] <- c.split$params

    tree[n.clust+2,] <- c(tree[id,1] + 1, 2*tree[id,2])

    Parent[n.clust+2] <- id

  }

  # if fewer than K clusters found, reset K

  if(length(ixs)<(2*K-1)){
    K <- (length(ixs)+1)/2
    split_indices <- split_indices[1:length(ixs)]
    tree <- tree[1:length(ixs),]
    vs <- vs[1:length(ixs),]
    bs <- bs[1:length(ixs)]
    dens <- dens[1:length(ixs)]
    rel.deps <- rel.deps[1:length(ixs)]
    Parent <- Parent[1:length(ixs)]
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
  for(i in 1:length(ixs)) Nodes[[i]] <- list(ixs = ixs[[i]], v = vs[i,], b = bs[i], params = pars[[i]], fval = dens[i], rel.dep = rel.deps[i], node = tree[i,], location = loci[i,])

  output <- list(cluster = asgn, model = tree, Parent = Parent, Nodes = Nodes, data = X, method = 'MDH', args = list(v0 = v0, minsize = minsize, bandwidth = bandwidth, split.index = split.index, alphamin = alphamin, alphamax = alphamax, maxit = maxit, ftol = ftol))

  class(output) <- 'ppci_cluster_solution'

  output
}


### function f_md() evaluates the projection index for MDPP. Assumes the data have zero mean
## arguments:
# v = projection vector
# X = data matrix
# P = list of parameters including (at least)
#   $h (bandwidth)
#   $alpha (constraint width. larger values allow hyperplanes further from the mean of the data)
#   $C (constant affecting the slope of the penalty)

## output is a scalar, the value of the density on the (best) hyperplane orthogonal to v

f_md <- function(v, X, P){

  # find density of data projected along v/||v||

  p <- X%*%v/norm_vec(v)
  n <- length(p)

  # if alpha is zero then hyperplane passes through origin

  if(P$alpha==0){
    return(sum(1/sqrt(2*pi)*exp(-p^2/2/P$h^2)/n/P$h))
  }

  # otherwise find the minimum density hyperplane along v within [-alpha*sigma, alpha*sigma]

  s <- sd(p)
  den <- density(p, bw = P$h, from = -P$alpha*s-1, to = P$alpha*s+1, n = 100)
  pen <- den$y + P$C*((den$x<(-P$alpha*s))*(-P$alpha*s-den$x)^2 + (den$x>(P$alpha*s))*(den$x-P$alpha*s)^2)
  mins <- which(apply(rbind(pen[1:98], pen[2:99], pen[3:100]), 2, function(x) x[2]<=min(x)))+1
  bs <- den$x[mins]
  fs <- pen[mins]

  # use bisection to refine the location of all local minima

  for(i in 1:length(mins)){
    hi <- den$x[mins[i]+1]
    lo <- den$x[mins[i]-1]
    repeat{
      mid <- (hi+lo)/2
      dmid <- dkde(mid, p, P$h, P$alpha, P$C, s)
      if(abs(dmid)<1e-8 || (hi-mid)<1e-10){
        bs[i] <- mid
        fs[i] <- 1/sqrt(2*pi)/n/P$h*sum(exp(-(p-mid)^2/2/P$h^2))
        if(bs[i]>(P$alpha*s)) fs[i] <- fs[i] + P$C*(bs[i]-P$alpha*s)^2
        if(bs[i]<(-P$alpha*s)) fs[i] <- fs[i] + P$C*(-bs[i]-P$alpha*s)^2
        break
      }
      else if(dmid<0) lo <- mid
      else hi <- mid
    }
  }

  # return minimum of the minima after refinement

  min(fs)
}

### function df_md() evaluates the gradient of the projection index for MDPP. Assumes data have zero mean
## arguments:
# v = projection vector
# X = data matrix
# P = list of parameters including (at least)
#   $h (bandwidth)
#   $alpha (constraint width. larger values allow hyperplanes further from the mean of the data)
#   $C (constant affecting the slope of the penalty)
#   $COV (covariance matrix of data)

## output is a vector. the gradient of the projection index at v

df_md = function(v, X, P){

  # first steps are as in evaluation of the projection index, to find the location of the optimum

  p <- X%*%v/norm_vec(v)
  n <- length(p)

  nv <- norm_vec(v)

  if(P$alpha==0){
    c(t(1/sqrt(2*pi)/n/P$h^3*(exp(-p^2/2/P$h^2)*(-p)))%*%(X/nv-(X%*%v)%*%t(v)/nv^3))
  }

  s <- sd(p)
  den <- density(p, bw = P$h, from = -P$alpha*s-1, to = P$alpha*s+1, n = 100)
  pen <- den$y + P$C*((den$x<(-P$alpha*s))*(-P$alpha*s-den$x)^2 + (den$x>(P$alpha*s))*(den$x-P$alpha*s)^2)
  mins <- which(apply(rbind(pen[1:98], pen[2:99], pen[3:100]), 2, function(x) x[2]<=min(x)))+1
  bs <- den$x[mins]
  fs <- pen[mins]
  for(i in 1:length(mins)){
    hi <- den$x[mins[i]+1]
    lo <- den$x[mins[i]-1]
    repeat{
      mid <- (hi+lo)/2
      dmid <- dkde(mid, p, P$h, P$alpha, P$C, s)
      if(abs(dmid)<1e-8 || (hi-mid)<1e-10){
        bs[i] <- mid
        fs[i] <- 1/sqrt(2*pi)/n/P$h*sum(exp(-(p-mid)^2/2/P$h^2))
        if(bs[i]>(P$alpha*s)) fs[i] <- fs[i] + P$C*(bs[i]-P$alpha*s)^2
        if(bs[i]<(-P$alpha*s)) fs[i] <- fs[i] + P$C*(-bs[i]-P$alpha*s)^2
        break
      }
      else if(dmid<0) lo <- mid
      else hi <- mid
    }
  }

  w <- which.min(fs)
  b <- bs[w]

  # compute the gradient of the the density on hyperplane H(v, b)

  if(b<(-P$alpha*s)){
    ds <- 1/nv^2*(P$COV%*%v/s-s*v)
    dpen <- c(-P$C*2*(-b-P$alpha*s)*P$alpha*ds)
    c(t(1/sqrt(2*pi)/n/P$h^3*(exp(-(p-b)^2/2/P$h^2)*(b-p)))%*%(X/nv-(X%*%v)%*%t(v)/nv^3)) + dpen
  }
  else if(b>(P$alpha*s)){
    ds <- 1/nv^2*(P$COV%*%v/s-s*v)
    dpen <- -c(P$C*2*(b-P$alpha*s)*P$alpha*ds)
    c(t(1/sqrt(2*pi)/n/P$h^3*(exp(-(p-b)^2/2/P$h^2)*(b-p)))%*%(X/nv-(X%*%v)%*%t(v)/nv^3)) + dpen
  }
  else{
    c(t(1/sqrt(2*pi)/n/P$h^3*(exp(-(p-b)^2/2/P$h^2)*(b-p)))%*%(X/nv-(X%*%v)%*%t(v)/nv^3))
  }
}

### function dkde() evaluates the gradient of the (penalised) kernel density estimator of the univariate
### projected dataset. Used in bisection to find the minimum density hyperplane along projection
## arguments:
# pt = point at which to evaluate the gradient
# xs = projected data points
# bw = bandwidth
# al = alpha parameter determining width of constraint
# K = constant affecting the slope of the penalty
# s = standard deviation of xs

## output is a scalar, the gradient of the projected density at pt

dkde <- function(pt, xs, bw, al, K, s){
  if(pt<(-al*s)) 1/sqrt(2*pi)/length(xs)/bw^3*sum(exp(-(xs-pt)^2/2/bw^2)*(xs-pt)) - 2*K*(-al*s-pt)
  else if(pt>(al*s)) 1/sqrt(2*pi)/length(xs)/bw^3*sum(exp(-(xs-pt)^2/2/bw^2)*(xs-pt)) + 2*K*(pt-al*s)
  else 1/sqrt(2*pi)/length(xs)/bw^3*sum(exp(-(xs-pt)^2/2/bw^2)*(xs-pt))
}

### function is_minim() checks if the optimal hyperplane along v is at a local minimum of the density
## arguments:
# v = projection vector
# X = data matrix
# P = list of parameters including (at least)
#   $h (bandwidth)
#   $alpha (constraint width. larger values allow hyperplanes further from the mean of the data)

## output is logical, is the best hyperplane orthogonal to v at a minimum of the density

is_minim <- function(v, X, P){

  # compute the projected density along v/||v||

  p <- X%*%v/norm_vec(v)
  s <- sd(p)
  d <- density(p, bw = P$h, n = 200)

  # if the minimum lies at one of the boundaries then it is not a local
  # minimum of the density

  pen <- d$y + P$C*((d$x<(-P$alpha*s))*(-P$alpha*s-d$x)^2 + (d$x>(P$alpha*s))*(d$x-P$alpha*s)^2)
  
  modes <- which(apply(rbind(d$y[3:200], d$y[2:199], d$y[1:198]), 2, function(x) x[2]>=max(x))) + 1
  bix <- which.min(pen)

  (bix > modes[1] && bix < max(modes) && min(sum(p<d$x[bix]), sum(p>d$x[bix]))>P$nmin)
}


### function mdpp() performs projection pursuit for finding minimal density hyperplanes
## arguments:
# v = initial projection vector
# X = data matrix
# P = list of parameters including (at least)
#   $h (bandwidth)
#   $alpha (constraint width. larger values allow hyperplanes further from the mean of the data)
#   $C (constant affecting the slope of the penalty)
#   $COV (covariance matrix of data)
# alphamin = initial constraint on the distance of hyperplane to the mean of the data
# alphamax = maximum allowable distance of hyperplane from the mean of the data
# verb = verbosity level. For values greater than 0 plots are produced to illustrate the progress of the algorithm
# labels = vector of class labels. used only in the plotting of progress
# maxit = maximum number of iterations in BFGS for each value of alpha
# ftol = tolerance for termination of BFGS based on function value changes

## output is a list containing the optimal projection vector and the corresponding value of alpha

mdpp <- function(v, X, P, alphamin, alphamax, verb, labels, maxit, ftol){

  # initialise with alpha = 0

  P$alpha <- alphamin
  v_opt <- v
  al_opt <- alphamin

  # increase alpha and store the solution for the largest value of alpha which is a local minimum

  while(P$alpha <= alphamax){
    v <- ppclust.optim(v, f_md, df_md, X, P, md_b, verbosity = verb, labels = labels, method = 'MDH', maxit = maxit, ftol = ftol)$par
    if(is_minim(v, X, P)){
      v_opt <- v
      al_opt <- P$alpha
    }
    P$alpha <- P$alpha + 0.1
  }
  list(v = v_opt/norm_vec(v_opt), alpha = al_opt)
}

### function md_b() determines the location of the optimal hyperplane along v
## arguments:
# v = projection vector
# X = data matrix
# P = list of parameters including (at least)
#   $h (bandwidth)
#   $alpha (constraint width. larger values allow hyperplanes further from the mean of the data)

## output is a scalar, the value of b making H(v, b) the minimum density hyperlpane orthogonal to v

md_b <- function(v, X, P){
  p <- X%*%v/norm_vec(v)
  n <- length(p)
  s <- sd(p)
  den <- density(p, bw = P$h, from = -P$alpha*s, to = P$alpha*s, n = 1000)
  w <- which.min(den$y)
  den$x[w[ceiling(length(w)/2)]]
}

### function mdh() determines the minimum density hyperplane
## arguments:
# X = dataset (matrix). each row is a datum. required
# v0 = initial projection direction(s). can be a matrix
#       in which each column is an initialisation to try.
#       can be a function of the data matrix (or subset
#       thereof corresponding to the cluster being split) which returns
#       a matrix in which each column is an initialisation.
#       optional, default is the first principal component
# bandwidth = positive numeric bandwidth parameter used in kernel density estimation (h).
#       optional, default is bandwidth = 0.9*eigen(cov(X))$values[1]^.5/nrow(X)^.2
# alphamin = initial width of constraint F(v) = [mu-alpha*sd, mu+alpha*sd]
# alphamax = maximum width of constraint F(v) = [mu - alpha*sd, mu+alpha*sd]
#       optional, default is 0.9
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

## output is a list of lists, the i-th stores the details of the minimum density hyperplane
## arising from the initialisation at v0[,i]. Each mdh has contains
# $cluster = the cluster assignment vector
# $v = the optimal projection vector
# $b = the value of b making H(v, b) the mdh
# rel.dep = the relative depth of H(v, b)
# fval = the density on H(v, b)
# params = list of parameters used to find H(v, b)

mdh <- function(X, v0 = NULL, minsize = NULL, bandwidth = NULL, alphamin = NULL, alphamax = NULL, verb = NULL, labels = NULL, maxit = NULL, ftol = NULL){

  if(is.data.frame(X)) X <- as.matrix(X)

  params <- list()

  if(is.null(minsize)) params$nmin <- 1
  else if(is.numeric(minsize)) params$nmin <- minsize
  else stop('minsize must be integer.')

  # if the data contain fewer than 2*minsize points, do not split

  if(is.vector(X)) n <- 1
  else n <- nrow(X)
  if(n<(2*params$nmin)){
    v <- numeric(ncol(X)) + 1/sqrt(ncol(X))
    b <- mean(X%*%v)
    rel.dep <- 0
    fval <- Inf
    cluster <- numeric(nrow(X)) + 1
    fitted <- X[,1:2]
    return(list(cluster = cluster, v = v, b = b, rel.dep = 0, fval = Inf, params = list(alpha=0,K=1000,h=1), fitted = fitted, data = X, method = 'MDH'))
  }



  # if labels are supplied, ensure they are integers 1:K (K the number of classes)

  if(!is.null(labels)){
    lab_new <- numeric(length(labels))
    u <- unique(labels)
    for(i in 1:length(u)) lab_new[which(labels==u[i])] = i
    labels <- lab_new
  }

  if(is.null(verb)) verb <- 0

  if(is.null(maxit)) maxit <- 50

  if(is.null(ftol)) ftol <- 1e-8

  # if the data do not have zero mean, center them

  mns <- colMeans(X)
  if(max(abs(mns))>1e-7) X <- t(t(X)-mns)

  # set up parameters for optimisation

  if(is.null(alphamin)) alphamin <- 0

  if(is.null(alphamax)) alphamax <- 1

  params$COV <- cov(X)

  if(is.null(bandwidth)){
    if(ncol(X)>2) params$h <- 0.9*rARPACK::eigs_sym(params$COV, 1)$values[1]^.5/n^0.2
    else params$h <- 0.9*eigen(params$COV)$values[1]^.5/n^0.2
  }
  else if(is.numeric(bandwidth)) params$h <- bandwidth
  else if(is.function(bandwidth)) params$h <- bandwidth(X)
  else stop('bandwidth must be numeric or a function of the data being split')

  params$C <- 100*exp(-0.5)/sqrt(2*pi)/params$h^2

  if(is.null(v0)){
    if(ncol(X)>2) E <- matrix(rARPACK::eigs_sym(params$COV, 1)$vectors, ncol = 1)
    else E <- matrix(eigen(params$COV)$vectors[,1], ncol = 1)
  }
  else if(is.function(v0)) E <- v0(X)
  else if(is.vector(v0)) E <- matrix(v0, ncol = 1)
  else E <- v0

  hyperplanes <- list()

  # find the mdh arising from each column of E (v0)

  for(i in 1:ncol(E)){

    v <- mdpp(E[,i], X, params, alphamin, alphamax, verb, labels, maxit, ftol)

    params$alpha <- v$alpha

    v <- v$v

    b <- md_b(v, X, params)

    pass <- X%*%v<b

    if(is_minim(v, X, params)){
      fval <- f_md(v, X, params)
      depth <- md_reldepth(v, X, params)
    }
    else{
      fval <- Inf
      depth <- 0
    }

    if(ncol(X)>2) v2 <- rARPACK::eigs_sym(cov(X-X%*%v%*%t(v)), 1)$vectors
    else v2 <- eigen(cov(X-X%*%v%*%t(v)))$vectors[,1]

    hyperplanes[[i]] <- list(cluster = pass+1, v = v, b = b + (mns%*%v)[1], rel.dep = depth, fval = fval, params = params, method = 'MDH', data = t(t(X)+mns), fitted = t(t(X)+mns)%*%cbind(v, v2))

    class(hyperplanes[[i]]) <- 'ppci_hyperplane_solution'
  }

  best_sol <- which.max(unlist(lapply(hyperplanes, function(l) l$rel.dep)))

  output <- hyperplanes[[best_sol]]

  output$alternatives <- hyperplanes[-best_sol]

  output
}

### function md_reldepth() computes the relative depth of the best hyperplane orthogonal to v
## arguments:
# v = projection vector
# X = data matrix
# P = list of parameters including (at least)
#   $h (bandwidth)
#   $alpha (constraint width. larger values allow hyperplanes further from the mean of the data)

## output is a scalar, the relative depth of the hypeprlane

md_reldepth <- function(v, X, P){

  # compute the projected density on v/||v|| and find the minimiser

  p <- X%*%v/norm_vec(v)
  n <- length(p)
  s <- sd(p)
  den <- density(p, bw = P$h, from = -P$alpha*s, to = P$alpha*s, n = 100)

  xmin <- den$x[which.min(den$y)]

  denmin <- min(den$y)+1e-10

  # compute the density to the left and right of the minimiser and find the relative depth

  den.left <- density(p, bw = P$h, to = xmin, n = 100)

  den.right <- density(p, bw = P$h, from = xmin, n = 100)

  (min(max(den.left$y), max(den.right$y))-denmin)/denmin
}



