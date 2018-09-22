mddr <- function(X, p, minsize = NULL, v0 = NULL, bandwidth = NULL, alphamin = NULL, alphamax = NULL, verb = NULL, labels = NULL, maxit = NULL, ftol = NULL){

  if(is.data.frame(X)) X <- as.matrix(X)

  if(is.null(verb)) verb = 0
  # set parameters for projection pursuit

  n <- nrow(X)
  d <- ncol(X)

  # stores hyperplane separators for each node, v and b

  vs <- matrix(0, d, p)

  # repeatedly identify the next projection dimension, based on orthogonalising the data to the projections already found

  for(it in 1:p){
    if(it==1) proj <- mdh(X, v0, minsize, bandwidth, alphamin, alphamax, verb, labels, maxit, ftol)
    else proj <- mdh(X-X%*%vs[,1:(it-1)]%*%t(vs[,1:(it-1)]), v0, minsize, bandwidth, alphamin, alphamax, verb, labels, maxit, ftol)
    vs[,it] <- proj$v
  }

  output <- list(projection = vs, fitted = X%*%vs, data = X, method = 'MDH', args = list(v0 = v0, minsize = minsize, bandwidth = bandwidth, alphamin = alphamin, alphamax = alphamax, maxit = maxit, ftol = ftol))

  class(output) <- 'ppci_projection_solution'

  output
}


mcdr <- function(X, p, minsize = NULL, v0 = NULL, verb = NULL, labels = NULL, maxit = NULL, ftol = NULL){

  if(is.data.frame(X)) X <- as.matrix(X)

  if(is.null(verb)) verb = 0
  # set parameters for projection pursuit

  n <- nrow(X)
  d <- ncol(X)

  # stores hyperplane separators for each node, v and b

  vs <- matrix(0, d, p)

  # repeatedly identify the next projection dimension, based on orthogonalising the data to the projections already found

  for(it in 1:p){
    if(it==1) proj <- mch(X, v0, minsize, verb, labels, maxit, ftol)
    else proj <- mch(X-X%*%vs[,1:(it-1)]%*%t(vs[,1:(it-1)]), v0, minsize, verb, labels, maxit, ftol)
    vs[,it] <- proj$v
  }

  output <- list(projection = vs, fitted = X%*%vs, data = X, method = 'MCDC', args = list(v0 = v0, minsize = minsize, maxit = maxit, ftol = ftol))

  class(output) <- 'ppci_projection_solution'

  output
}



ncutdr <- function(X, p, v0 = NULL, s = NULL, minsize = NULL, verb = NULL, labels = NULL, maxit = NULL, ftol = NULL){

  if(is.data.frame(X)) X <- as.matrix(X)

  if(is.null(verb)) verb = 0
  # set parameters for projection pursuit

  n <- nrow(X)
  d <- ncol(X)

  # stores hyperplane separators for each node, v and b

  vs <- matrix(0, d, p)

  # repeatedly identify the next projection dimension, based on orthogonalising the data to the projections already found

  for(it in 1:p){
    if(it==1) proj <- ncuth(X, v0, s, minsize, verb, labels, maxit, ftol)
    else proj <- ncuth(X-X%*%vs[,1:(it-1)]%*%t(vs[,1:(it-1)]), v0, s, minsize, verb, labels, maxit, ftol)
    vs[,it] <- proj$v
  }

  output <- list(projection = vs, fitted = X%*%vs, data = X, method = 'NCutH', args = list(v0 = v0, s=s, minsize = minsize, maxit = maxit, ftol = ftol))

  class(output) <- 'ppci_projection_solution'

  output
}
