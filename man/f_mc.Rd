\name{f_mc}
\alias{f_mc}
\title{Variance Ratio Clusterability Across a Hyperplane}
\description{
  Finds the variance ratio clusterability measured across the best hyperplane orthogonal to a given projection vector. Used as projection index for projection pursuit to find maximum clusterability hyperplanes.
}
\usage{
f_mc(v, X, P)
}
\arguments{
  \item{v}{a numeric vector of length ncol(X)}
  \item{X}{a numeric matrix (num_data x num_dimensions) to be projected on v}
  \item{P}{a list of parameters including (at least) $nmin (positive integer minimum cluster size)}
}
\value{
  the value of the variance ratio across the best hyperplane orthogonal to v
}
\keyword{internal}
