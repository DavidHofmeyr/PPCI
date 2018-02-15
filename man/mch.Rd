\name{mch}
\alias{mch}
\title{Maximum Clusteriability Hyperplane}
\description{
  Finds maximum clusterability hyperplane(s) for clustering.
}
\usage{
mch(X, v0, minsize, verb, labels, maxit, ftol)
}
\arguments{
  \item{X}{a numeric matrix (num_data x num_dimensions); the dataset to be clustered.}
  \item{v0}{(optional) initial projection direction(s). a matrix with ncol(X) rows. each column of v0 is used as an initialisation for projection pursuit. if omitted then a single initialisation is used; the vector joining the cluster means from a 2-means solution.}
  \item{minsize}{(optional) the minimum cluster size allowable. if omitted then minsize = 1.}
  \item{verb}{(optional) verbosity level of optimisation procedure. verb==0 produces no output. verb==1 produces plots illustrating the progress of projection pursuit via plots of the projected data. verb==2 adds to these plots additional information about the progress. verb==3 creates a folder in working directory and stores all plots for verb==2. if omitted then verb==0.}
  \item{labels}{(optional) vector of class labels. not used in the actual clustering procedure. only used for illustrative purposes for values of verb>0.}
  \item{maxit}{(optional) maximum number of iterations in optimisation. if omitted then maxit=50.}
  \item{ftol}{(optional) tolerance level for convergence of optimisation, based on relative function value improvements. if omitted then ftol = 1e-8.}
}
\value{
  an unnamed list containing one entry for each initialisation (column of v0). each entry is a named list with the following components
  \item{$cluster}{cluster assignment vector.}
  \item{$v}{the optimal projection vector.}
  \item{$b}{the value of b making H(v, b) the minimum normalised cut hyperplane.}
  \item{$fval}{the variance ratio clusterability across H(v, b).}
  \item{$method}{=="MCDC".}
  \item{$params}{list of parameters used to find H(v, b).}
}
\references{
  Hofmeyr, D., Pavlidis, N. (2015) Maximum Clusterability Divisive Clustering. \emph{Computational Intelligence, 2015 IEEE Symposium Series on}, pp. 780--786.
}
\examples{
## generate dataset with elongated clusters for which variance ratio in
## both dimensions is misleading for clustering
set.seed(1)
S <- matrix(c(1, .7, .7, 1), 2, 2)
X <- matrix(rnorm(2000), ncol = 2)\%*\%S
X[,1] <- X[,1] + rep(c(.8, -.8), each = 500)
X[,2] <- X[,2] + rep(c(-.8, .8), each = 500)

## find the optimal variance ratio hyperplane solution
sol <- mch(X)

## visualise the solution
plot(X, col = sol[[1]]$cluster)
}
\keyword{file}
