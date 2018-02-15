\name{df_mc}
\alias{df_mc}
\title{Gradient of the Variance Ratio Clusterability Across a Hyperplane}
\description{
  Finds the gradient of the variance ratio across the best hyperplane orthogonal to a given projection vector. Used to obtain maximum clusterability hyperplanes using gradient based optimisation.
}
\usage{
df_mc(v, X, P)
}
\arguments{
  \item{v}{a numeric vector of length ncol(X)}
  \item{X}{a numeric matrix (num_data x num_dimensions) to be projected on v}
  \item{P}{a list of parameters including (at least) $nmin (positive integer minimum cluster size).}
}
\value{
  the (vector) gradient of the variance across the best hyperplane orthogonal to v.
}
\keyword{internal}
