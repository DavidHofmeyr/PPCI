\name{mc_b}
\alias{mc_b}
\title{Location of Optimal Variance Ratio Hyperplane}
\description{
  The value of b which maximises the variance ratio across the hyperplane H(v, b) = [x : v'x = b]. If v is a locally optimal solution for the projection index f_mc then H(v, b) is a locally optimal hyperplane.
}
\usage{
mc_b(v, X, P)
}
\arguments{
  \item{v}{a numeric vector of length ncol(X)}
  \item{X}{a numeric matrix (num_data x num_dimensions) to be projected on v}
  \item{P}{a list of parameters including (at least) $nmin (positive integer minimum cluster size).}
}
\value{
  the value of b given in the description.
}
\keyword{internal}
