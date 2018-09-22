\name{ncutdc}
\alias{ncutdc}
\title{Divisive Clustering Using Minimum Normalised Cut Hyperplanes}
\description{
  Generates a binary partitioning tree by recursively partitioning a dataset using a hierarchical collection of hyperplanes with
	low normalised cut measured across them.
}
\usage{
ncutdc(X, K, split.index, v0, s, minsize, verb, labels, maxit, ftol)
}
\arguments{
  \item{X}{a numeric matrix (num_data x num_dimensions); the dataset to be clustered.}
  \item{K}{the number of clusters to extract.}
  \item{split.index}{(optional) determines the order in which clusters are split (in decreasing order of split indices). can be a numeric valued function(v, X, P) of projection vector v, data matrix X and list of parameters P. can also be one of "size" (split the largest cluster) or "fval" (split the cluster with the minimum density hyperplane). if omitted then "fval" is used.}
  \item{v0}{(optional) initial projection direction(s). a function(X) of the data being split, which returns a matrix with ncol(X) rows. each column of the output of v0(X) is used as an initialisation for projection pursuit. the solution with the minimum normalised cut is used within the final model. initialisations are determined separately for each cluster being split. if omitted then a single initialisation is used; the first principal component.}
  \item{s}{(optional) used to compute the scaling parameter (sigma) for pairwise similarities. a numeric valued function(X) of the cluster being split. if omitted then s(X) = 100*eigen(cov(X))$values[1]^.5*nrow(X)^(-0.2).}
  \item{minsize}{(optional) the minimum cluster size allowable. if omitted then minsize = 1.}
  \item{verb}{(optional) verbosity level of optimisation procedure. verb==0 produces no output. verb==1 produces plots illustrating the progress of projection pursuit via plots of the projected data. verb==2 adds to these plots additional information about the progress. verb==3 creates a folder in working directory and stores all plots for verb==2. if omitted then verb==0.}
  \item{labels}{(optional) vector of class labels. not used in the actual clustering procedure. only used for illustrative purposes for values of verb>0.}
  \item{maxit}{(optional) maximum number of iterations in optimisation. if omitted then maxit=50.}
  \item{ftol}{(optional) tolerance level for convergence of optimisation, based on relative function value improvements. if omitted then ftol = 1e-8.}
}
\value{
  a named list with class ppci_cluster_solution containing
  \item{$cluster}{cluster assignment vector.}
  \item{$model}{matrix containing the would-be location of each node (depth and position at depth) within a complete tree of appropriate depth.}
  \item{$nodes}{unnamed list each element of which is a named list containing details of the binary partitions at each node in the model.}
  \item{$data}{the data matrix being clustered.}
  \item{$method}{=="NCutH". used in plotting and model modification functions.}
  \item{$args}{named list of arguments passed to ncutdc.}
}
\references{
  Hofmeyr, D. (2016) Clustering by Minimum Cut Hyperplanes. \emph{IEEE Transactions on Pattern Analysis and Machine Intelligence}, \bold{39}(8), 1547 -- 1560.
}
\examples{
## load dermatology dataset
data(dermatology)

## obtain clustering solution
sol <- ncutdc(dermatology$x, 6)

## evaluate solution using external cluster validity metrics
cluster_performance(sol$cluster, dermatology$c)
}
\keyword{file}
