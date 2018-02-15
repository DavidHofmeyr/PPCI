\name{df_md}
\alias{df_md}
\title{Gradient of the Integrated Density on a Hyperplane}
\description{
  Finds the gradient of the integrated density of the best hyperplanes orthogonal to a given projection vector (assumes the data have zero mean vector). Used to obtain minimum density hyperplanes using gradient based optimisation.
}
\usage{
df_md(v, X, P)
}
\arguments{
  \item{v}{a numeric vector of length ncol(X)}
  \item{X}{a numeric matrix (num_data x num_dimensions) to be projected on v}
  \item{P}{a list of parameters including (at least) $h (positive numeric bandwidth value), $alpha (positive numeric constraint width), $C (positive numeric affecting the slope of the penalty), $COV (covariance matrix of X)}
}
\value{
  the (vector) gradient of the integrated density of the best hyperplane orthogonal to v.
}
\keyword{internal}
