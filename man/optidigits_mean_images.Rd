\name{optidigits_mean_images}
\alias{optidigits_mean_images}
\title{Visualise Cluster Means from optidigits data set}
\description{
  Provides a visualisation of the cluster means computed from the optidigits data set, recast as images. Cluster labels are aligned with the true labels using simulated annealing to maximise the trace of the confusion matrix (or subset if
  number of clusters != number of classes.)
}
\usage{
optidigits_mean_images(clusters)
}
\arguments{
  \item{clusters}{a vector of cluster assignments. Must take values in 1:k, where k is the number of clusters.}
}
\references{
  Lichman, M. (2013) UCI Machine Learning Repository. Irvine, CA: University of California, School of Information and Computer Science. \url{https://archive.ics.uci.edu/ml}
}
\examples{
### not run
run = FALSE
if(run){
  ## load optidigits dataset
  data(optidigits)

  ## obtain a clustering solution using normalised cut hyperplanes
  sol <- ncutdc(optidigits$x, 10)

  ## visualise the cluster means as images
  optidigits_mean_images(sol$cluster)
  }
}
\keyword{file}
