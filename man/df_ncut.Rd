\name{df_ncut}
\alias{df_ncut}
\title{Gradient of the Normalised Cut Across a Hyperplane}
\description{
  Finds the gradient of the normalised cut across the best hyperplane orthogonal to a given projection vector. Used to obtain minimum normalised cut hyperplanes using gradient based optimisation.
}
\usage{
df_ncut(v, X, P)
}
\arguments{
  \item{v}{a numeric vector of length ncol(X)}
  \item{X}{a numeric matrix (num_data x num_dimensions) to be projected on v}
  \item{P}{a list of parameters including (at least) $s (positive numeric scaling parameter), $nmin (positive integer minimum cluster size).}
}
\value{
  the (vector) gradient of the normalised cut across the best hyperplane orthogonal to v.
}
\keyword{internal}
