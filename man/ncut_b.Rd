\name{ncut_b}
\alias{ncut_b}
\title{Location of Minimum Normalised Cut Hyperplane}
\description{
  The value of b which minimises the normalised cut across the hyperplane H(v, b) = [x : v'x = b]. If v is a locally optimal solution for the projection index f_ncut then H(v, b) is a locally optimal hyperplane.
}
\usage{
ncut_b(v, X, P)
}
\arguments{
  \item{v}{a numeric vector of length ncol(X)}
  \item{X}{a numeric matrix (num_data x num_dimensions) to be projected on v}
  \item{P}{a list of parameters including (at least) $s (positive numeric scaling parameter), $nmin (positive integer minimum cluster size).}
}
\value{
  the value of b given in the description.
}
\keyword{internal}
