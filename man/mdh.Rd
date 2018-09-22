\name{mdh}
\alias{mdh}
\title{Minimum Density Hyperplane}
\description{
  Finds minimum density hyperplane(s) for clustering.
}
\usage{
mdh(X, v0, minsize, bandwidth, alphamin, alphamax, verb, labels, maxit, ftol)
}
\arguments{
  \item{X}{a numeric matrix (num_data x num_dimensions); the dataset to be clustered.}
  \item{v0}{(optional) initial projection direction(s). a matrix with ncol(X) rows. each column of v0 is used as an initialisation for projection pursuit. if omitted then a single initialisation is used; the first principal component.}
  \item{minsize}{(optional) minimum cluster size. if omitted then minsize = 1.}
  \item{bandwidth}{(optional) positive numeric bandwidth parameter (h) for the kernel density estimator. if omitted then bandwidth = 0.9*eigen(cov(X))$values[1]^.5*nrow(X)^(-0.2).}
  \item{alphamin}{(optional) initial (scaled) bound on the distance of the optimal hyperplane from the mean of the data. if omitted then alphamin = 0.}
  \item{alphamax}{(optional) maximum/final (scaled) distance of the optimal hyperplane from the mean of the data. if omitted then alphamax = 1.}
  \item{verb}{(optional) verbosity level of optimisation procedure. verb==0 produces no output. verb==1 produces plots illustrating the progress of projection pursuit via plots of the projected data. verb==2 adds to these plots additional information about the progress. verb==3 creates a folder in working directory and stores all plots for verb==2. if omitted then verb==0.}
  \item{labels}{(optional) vector of class labels. not used in the actual clustering procedure. only used for illustrative purposes for values of verb>0.}
  \item{maxit}{(optional) maximum number of iterations in optimisation for each value of alpha. if omitted then maxit=15.}
  \item{ftol}{(optional) tolerance level for convergence of optimisation, based on relative function value improvements. if omitted then ftol = 1e-5.}
}
\value{
  a named list with class ppci_hyperplane_solution with the following components
  \item{$cluster}{cluster assignment vector.}
  \item{$v}{the optimal projection vector.}
  \item{$b}{the value of b making H(v, b) the minimum density hyperplane.}
  \item{$fitted}{data projected into two dimensional subspace defined by $v and the principal component in the null space of $v.}
  \item{$data}{the input data matrix.}
  \item{$rel.dep}{the relative depth of H(v, b).}
  \item{$fval}{the integrated dentsity on H(v, b).}
  \item{$method}{=="MDH".}
  \item{$params}{list of parameters used to find H(v, b).}
  \item{$alternatives}{an unnamed list. If more than one initilisation is considered, the alternatives to the best are stored in this field.}
}
\references{
  Pavlidis N.G., Hofmeyr D.P., Tasoulis S.K. (2016) Minimum Density Hyperplanes. \emph{Journal of
Machine Learning Research}, \bold{17}(156), 1--33.
}
\examples{
## load breast cancer dataset
data(breastcancer)

## find minimum density hyperplane
sol <- mdh(breastcancer$x)

## visualise the solution
plot(sol)

## evaluate the quality of the partition
success_ratio(sol$cluster, breastcancer$c)
}
\keyword{file}
