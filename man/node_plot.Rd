\name{node_plot}
\alias{node_plot}
\title{Visualise a Node in a Hierarchical Clustering Model}
\description{
  Provides a visualisation of the partition at an internal node (or of the potential partition at a leaf node) via a two-dimensional projection of the data assigned to the node.
}
\usage{
node_plot(sol, node, labels)
}
\arguments{
  \item{sol}{a clustering solution arising from one of the functions mcdc, mddc and ncutdc.}
  \item{node}{the node to be visualised. can be either an integer specifying the node number in sol$Nodes or a vector of length two specifying c(depth, position at depth) of the node.}
  \item{labels}{(optional) a vector of class labels. if provided then points in different classes are plotted in different colours. the external validity of each partition is also indicated when this can be computed using true class labels.}
}
\examples{
## load optidigits dataset
data(optidigits)

## generate a divisive clustering model using NCUTDC
sol <- ncutdc(optidigits$x, 10)

## visualise the node at depth 3 and position 2
node_plot(sol, c(3, 3))
}
\keyword{file}
