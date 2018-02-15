\name{f_ncut}
\alias{f_ncut}
\title{Normalised Cut Across a Hyperplane}
\description{
  Finds the normalised cut across the best hyperplane orthogonal to a given projection vector. Used as projection index for projection pursuit to find minimum normalised cut hyperplanes.
}
\usage{
f_ncut(v, X, P)
}
\arguments{
  \item{v}{a numeric vector of length ncol(X)}
  \item{X}{a numeric matrix (num_data x num_dimensions) to be projected on v}
  \item{P}{a list of parameters including (at least) $s (positive numeric scaling parameter), $nmin (positive integer minimum cluster size)}
}
\value{
  the value of the normalised cut across the best hyperplane orthogonal to v
}
\keyword{internal}
