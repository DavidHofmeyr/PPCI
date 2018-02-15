\name{f_md}
\alias{f_md}
\title{Integrated Density on a Hyperplane}
\description{
  Finds the integrated density of the best hyperplanes orthogonal to a given projection vector (assumes the data have zero mean vector). Used as projection index for projection pursuit to find minimum density hyperplanes.
}
\usage{
f_md(v, X, P)
}
\arguments{
  \item{v}{a numeric vector of length ncol(X)}
  \item{X}{a numeric matrix (num_data x num_dimensions) to be projected on v}
  \item{P}{a list of parameters including (at least) $h (positive numeric bandwidth value), $alpha (positive numeric constraint width), $C (positive numeric affecting the slope of the penalty)}
}
\value{
  the value of the integrated density of the best hyperplane orthogonal to v
}
\keyword{internal}
