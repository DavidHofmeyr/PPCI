\name{mddr}
\alias{mddr}
\title{Minimum Density Dimension Reduction}
\description{
  Finds a linear projection of a data set using projection pursuit to find the vector(s) orthogonal to minimum density hyperplanes.
}
\usage{
mddr(X, p, minsize, v0, bandwidth,
      alphamin, alphamax, verb, labels, maxit, ftol)
}
\arguments{
  \item{X}{a numeric matrix (num_data x num_dimensions); the dataset.}
  \item{p}{an integer; the number of dimensions in the projection.}
  \item{v0}{(optional) initial projection direction(s). a function(X) of the data, which returns a matrix with ncol(X) rows. each column of the output of v0(X) is used as an initialisation for projection pursuit. the solution with the minimum normalised cut is used within the final model. if omitted then a single initialisation is used for each column of the projection matrix; the first principal component within the null space of the other columns.}
  \item{bandwidth}{(optional) used to compute the bandwidth parameter (h) for the kernel density estimator. a numeric valued function(X) of the cluster being split. if omitted then bandwidth(X) = 0.9*eigen(cov(X))$values[1]^.5*nrow(X)^(-0.2).}
  \item{alphamin}{(optional) initial (scaled) bound on the distance of the optimal hyperplane from the mean of the data. if omitted then alphamin = 0.}
  \item{alphamax}{(optional) maximum (scaled) distance of the optimal hyperplane from the mean of the data. if omitted then alphamax = 1.}
  \item{verb}{(optional) verbosity level of optimisation procedure. verb==0 produces no output. verb==1 produces plots illustrating the progress of projection pursuit via plots of the projected data. verb==2 adds to these plots additional information about the progress. verb==3 creates a folder in working directory and stores all plots for verb==2. if omitted then verb==0.}
  \item{labels}{(optional) vector of class labels. not used in the actual projection pursuit. only used for illustrative purposes for values of verb>0.}
  \item{maxit}{(optional) maximum number of iterations in optimisation for each value of alpha. if omitted then maxit=15.}
  \item{ftol}{(optional) tolerance level for convergence of optimisation, based on relative function value improvements. if omitted then ftol = 1e-5.}
  \item{minsize}{(optional) the minimum number of data on each side of a hyperplane. if omitted then minsize = 1.}
}
\value{
  a named list with class ppci_projection_solution with the following components
  \item{$projection}{the num_dimensions x p projection matrix.}
  \item{$fitted}{the num_data x p projected data set.}
  \item{$data}{the input data matrix.}
  \item{$method}{=="MDH".}
  \item{$params}{list of parameters used to find $projection.}
}
\references{
  Pavlidis N.G., Hofmeyr D.P., Tasoulis S.K. (2016) Minimum Density Hyperplanes. \emph{Journal of
Machine Learning Research}, \bold{17}(156), 1--33.
}
\examples{
### not run
run = FALSE
if(run){
  ## load optidigits dataset
  data(optidigits)

  ## find nine dimensional projection (one fewer than
  ## the number of clusters, as is common in clustering)
  sol <- mddr(optidigits$x, 9)

  ## visualise the solution via the first 3 pairs of dimensions
  plot(sol, pairs = 3, labels = optidigits$c)

  ## compare with PCA projection
  pairs(optidigits$x\%*\%eigen(cov(optidigits$x))$vectors[,1:3], col = optidigits$c)
  }
}
\keyword{file}
