\name{md_reldepth}
\alias{md_reldepth}
\title{Relative Depth of a Hyperplane}
\description{
  Finds the relative depth, min((M1-m)/m, (M2-m)/m), where M1 is the maximum density below the optimal hyperplane orthogonal to v and M2 the maximum above, while m is the integrated density on the optimal hyperplane orthogonal to v.
}
\usage{
md_reldepth(v, X, P)
}
\arguments{
  \item{v}{a numeric vector of length ncol(X)}
  \item{X}{a numeric matrix (num_data x num_dimensions) to be projected on v}
  \item{P}{a list of parameters including (at least) $h (positive numeric bandwidth value), $alpha (positive numeric constraint width), $C (positive numeric affecting the slope of the penalty)}
}
\value{
  the value of the relative depth of the best hyperplane orthogonal to v.
}
\keyword{internal}
