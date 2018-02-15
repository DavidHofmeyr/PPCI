\name{tree_plot}
\alias{tree_plot}
\title{Visualise a Hierarchical Clustering Model}
\description{
  Provides a visualisation of a hierarchical clustering model via two-dimensional projections of the data assigned to each node.
}
\usage{
tree_plot(sol, labels, node.numbers)
}
\arguments{
  \item{sol}{a clustering solution arising from one of the functions mcdc, mddc and ncutdc.}
  \item{labels}{(optional) a vector of class labels. if provided then points in different classes are plotted in different colours.}
  \item{node.numbers}{(optional) logical. if TRUE then numbers are added to the plot to indicate the order in which nodes were added to the model. if omitted then node.numbers = TRUE.}
}
\examples{
## load optidigits dataset
data(optidigits)

## generate hierarchical clustering model using NCUTDC algorithm
sol <- ncutdc(optidigits$x, 10)

## visualise the solution
tree_plot(sol)
}
\keyword{file}
