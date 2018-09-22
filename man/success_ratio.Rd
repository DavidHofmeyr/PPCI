\name{success_ratio}
\alias{success_ratio}
\title{Evaluate External Valifity os a Binary Partition}
\description{
  Computes the success ratio of a binary partition by comparing the solution with true class labels.
}
\usage{
success_ratio(assigned, labels)
}
\arguments{
  \item{assigned}{a vector of cluster assignments made by a clustering algorithm.}
  \item{labels}{a vector of true class labels to be compared with assigned.}
}
\value{
  the success ratio of the cluster assignment solution.
}
\references{
  Hofmeyr, D. (2016) Clustering by Minimum Cut Hyperplanes. \emph{IEEE Transactions on Pattern Analysis and Machine Intelligence}.
}
\examples{
## load optidigits dataset
data(optidigits)

## generate a binary partition using minimum normalised cut hyperplane
sol <- ncuth(optidigits$x)

## evaluate using success ratio
success_ratio(sol$cluster, optidigits$c)
}
\keyword{file}
