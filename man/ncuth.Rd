\name{ncuth}
\alias{ncuth}
\title{Minimum Normalised Cut Hyperplane}
\description{
  Finds minimum normalised cut hyperplane(s) for clustering.
}
\usage{
ncuth(X, v0, s, minsize, verb, labels, maxit, ftol)
}
\arguments{
  \item{X}{a numeric matrix (num_data x num_dimensions); the dataset to be clustered.}
  \item{v0}{(optional) initial projection direction(s). a matrix with ncol(X) rows. each column of v0 is used as an initialisation for projection pursuit. if omitted then a single initialisation is used; the first principal component.}
  \item{s}{(optional) positive numeric scaling parameter (sigma). if omitted then s = 100*eigen(cov(X))$values[1]^.5*nrow(X)^(-0.2).}
  \item{minsize}{(optional) the minimum cluster size allowable. if omitted then minsize = 1.}
  \item{verb}{(optional) verbosity level of optimisation procedure. verb==0 produces no output. verb==1 produces plots illustrating the progress of projection pursuit via plots of the projected data. verb==2 adds to these plots additional information about the progress. verb==3 creates a folder in working directory and stores all plots for verb==2. if omitted then verb==0.}
  \item{labels}{(optional) vector of class labels. not used in the actual clustering procedure. only used for illustrative purposes for values of verb>0.}
  \item{maxit}{(optional) maximum number of iterations in optimisation. if omitted then maxit=50.}
  \item{ftol}{(optional) tolerance level for convergence of optimisation, based on relative function value improvements. if omitted then ftol = 1e-8.}
}
\value{
  a named list with class ppci_hyperplane_solution with the following components
  \item{$cluster}{cluster assignment vector.}
  \item{$v}{the optimal projection vector.}
  \item{$b}{the value of b making H(v, b) the minimum normalised cut hyperplane.}
  \item{$fitted}{data projected into two dimensional subspace defined by $v and the principal component in the null space of $v.}
  \item{$data}{the input data matrix.}
  \item{$fval}{the normalised cut across H(v, b).}
  \item{$method}{=="NCutH".}
  \item{$params}{list of parameters used to find H(v, b).}
  \item{$alternatives}{an unnamed list. If more than one initilisation is considered, the alternatives to the best are stored in this field.}
}
\references{
  Hofmeyr, D. (2016) Clustering by Minimum Cut Hyperplanes. \emph{IEEE Transactions on Pattern Analysis and Machine Intelligence}, \bold{39}(8), 1547 -- 1560.
}
\examples{
## load breast cancer dataset
data(breastcancer)

## find minimum normalised cut hyperplane
sol <- ncuth(breastcancer$x)

## visualise the solution
plot(sol)

## evaluate the performance of the solution
success_ratio(sol$cluster, breastcancer$c)
}
\keyword{file}
