\name{plot.ppci_cluster_solution}
\alias{plot.ppci_cluster_solution}
\title{Visualise a Hierarchical Clustering Model, or a Node Within a Hierarchical Clustering Model}
\description{
  Provides a visualisation of a hierarchical clustering model via two-dimensional projections of the data assigned to each node. Alternatively a more detailed visualisation of a single node in a hierarchical model.
}
\usage{
\method{plot}{ppci_cluster_solution}(x, node = NULL, labels = NULL, node.numbers = NULL, transparency = NULL, ...)
}
\arguments{
  \item{x}{a clustering solution arising from one of the functions mcdc, mddc and ncutdc.}
  \item{node}{(optional) integer corresponding to a node number (see node.numbers below). If provided then a detailed plot of the specified node is provided, instead of the entire cluster hierarchy.}
  \item{labels}{(optional) a vector of class labels. if provided then points in different classes are plotted in different colours.}
  \item{node.numbers}{(optional) logical. if TRUE then numbers are added to the plot to indicate the order in which nodes were added to the model. if omitted then node.numbers = TRUE.}
  \item{transparency}{(optional) if ommitted then points in scatterplots are
  shown as solid. If set to a value in (0, 1) then points are shown with transparency. Lower values correspond with a greater degree of transparency.}
  \item{...}{additional graphical parameters. Currently none implemented.}
}
\seealso{mcdc, mddc, ncutdc, tree_plot, node_plot}
\value{None}
\keyword{internal}
