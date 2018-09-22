\name{ncutdr}
\alias{ncutdr}
\title{Minimum Normalised Cut Dimension Reduction}
\description{
  Finds a linear projection of a data set using projection pursuit based on minimising the normalised cut measured in each dimension separately.
}
\usage{
ncutdr(X, p, v0, s, minsize, verb, labels, maxit, ftol)
}
\arguments{
  \item{X}{a numeric matrix (num_data x num_dimensions); the dataset.}
  \item{p}{an integer; the number of dimensions in the projection.}
  \item{v0}{(optional) initial projection direction(s). a function(X) of the data, which returns a matrix with ncol(X) rows. each column of the output of v0(X) is used as an initialisation for projection pursuit. the solution with the minimum normalised cut is used within the final model. if omitted then a single initialisation is used for each column of the projection matrix; the first principal component within the null space of the other columns.}
  \item{s}{(optional) used to compute the scaling parameter (sigma) for pairwise similarities. a numeric valued function(X) of the data. if omitted then s(X) = 100*eigen(cov(X))$values[1]^.5*nrow(X)^(-0.2).}
  \item{minsize}{(optional) the minimum cluster size allowable when computing the normalised cut. if omitted then minsize = 1.}
  \item{verb}{(optional) verbosity level of optimisation procedure. verb==0 produces no output. verb==1 produces plots illustrating the progress of projection pursuit via plots of the projected data. verb==2 adds to these plots additional information about the progress. verb==3 creates a folder in working directory and stores all plots for verb==2. if omitted then verb==0.}
  \item{labels}{(optional) vector of class labels. not used in the actual projection pursuit. only used for illustrative purposes for values of verb>0.}
  \item{maxit}{(optional) maximum number of iterations in optimisation. if omitted then maxit=50.}
  \item{ftol}{(optional) tolerance level for convergence of optimisation, based on relative function value improvements. if omitted then ftol = 1e-8.}
}
\value{
  a named list with class ppci_projection_solution with the following components
  \item{$projection}{the num_dimensions x p projection matrix.}
  \item{$fitted}{the num_data x p projected data set.}
  \item{$data}{the input data matrix.}
  \item{$method}{=="NCutH".}
  \item{$params}{list of parameters used to find $projection.}
}
\references{
  Hofmeyr, D. (2016) Clustering by Minimum Cut Hyperplanes. \emph{IEEE Transactions on Pattern Analysis and Machine Intelligence}, \bold{39}(8), 1547 -- 1560.
}
\examples{
### not run
run = FALSE
if(run){
  ## load optidigits dataset
  data(optidigits)

  ## find nine dimensional projection (one fewer than
  ## the number of clusters, as is common in clustering)
  sol <- ncutdr(optidigits$x, 9)

  ## visualise the solution via the first 3 pairs of dimensions
  plot(sol, pairs = 3, labels = optidigits$c)

  ## compare with PCA projection
  pairs(optidigits$x\%*\%eigen(cov(optidigits$x))$vectors[,1:3], col = optidigits$c)
  }
}
\keyword{file}
