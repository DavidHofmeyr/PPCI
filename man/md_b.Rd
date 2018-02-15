\name{md_b}
\alias{md_b}
\title{Location of Minimum Density Hyperplane}
\description{
  The value of b which minimises the penalised integrated density on the hyperplane H(v, b) = [x : v'x = b] (assumes the data have zero mean vector). If v is a locally optimal solution for the projection index f_md then H(v, b) is a locally optimal hyperplane.
}
\usage{
md_b(v, X, P)
}
\arguments{
  \item{v}{a numeric vector of length ncol(X)}
  \item{X}{a numeric matrix (num_data x num_dimensions) to be projected on v}
  \item{P}{a list of parameters including (at least) $h (positive numeric bandwidth value), $alpha (positive numeric constraint width), $C (positive numeric affecting the slope of the penalty)}
}
\value{
  the value of b given in the description.
}
\keyword{internal}
