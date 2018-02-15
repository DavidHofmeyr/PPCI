\name{tree_split}
\alias{tree_split}
\title{Split a Leaf in a Hierarchical Clustering Model}
\description{
  Adds an additional binary partition to an existing hierarchical clustering model produced by one of mcdc, mddc and ncutdc.
}
\usage{
tree_split(sol, node, ...)
}
\arguments{
  \item{sol}{a clustering solution arising from one of the functions mcdc, mddc and ncutdc.}
  \item{node}{the node to be further partitioned. can be either an integer specifying the node number in sol$nodes or a vector of length two specifying c(depth, position at depth) of the node.}
    \item{...}{any modifications to parameters used in optimisation. these should have the same names and types as the corresponding arguments for the method used to construct sol.}
}
\value{
  a list with the same components as sol. the $args field will reflect any changes included in ... above.
}
\examples{
## load the optidigits dataset
data(optidigits)

## cluster using minimum normalised cut hyperplanes,
## assuming no domain knowledge begin with 8 clusters
sol <- ncutdc(optidigits$x, 8)

## visualise solution
tree_plot(sol, node.numbers = TRUE)

## node 13 shows evidence of multiple clusters. Inspect this node more closely
node_plot(sol, 13)

## split node 13
sol_new <- tree_split(sol, 13)

## compare the solutions using external cluster validity metrics
cluster_performance(sol$cluster, optidigits$c)

cluster_performance(sol_new$cluster, optidigits$c)
}
\keyword{file}
