\name{PPCI-package}
\alias{PPCI-package}
\alias{PPCI}
\docType{package}
\title{
\packageTitle{PPCI}
}
\description{
This package provides implementations of three recently developed projection pursuit methods for clustering. These methods optimise measures of clusterability of the univariate (projected) dataset that are motivated by three well established
approaches to clustering; namely density clustering, centroid based clustering and clustering by graph cuts. The resulting partitions are formed by hyperplanes orthogonal to the optimal projection vectors, and multiple such partitioning hyperplanes are combined in a hierarchical model to generate complete clustering solutions. Model visualisations through low dimensional projections of the data/clusters are provided through multiple plotting functions, which facilitate model validation. Simple model modification functions then allow for pseudo-interactive clustering.

The three main clustering algorithms are implemented in the functions mddc, mcdc and ncutdc. Each takes two mandatory arguments, the data matrix (X) and the number of clusters (K). Numerous optional arguments allow the user to modify the specifics of optimisation, etc. If the correct number of clusters is not known an approximate number can be used, and the resulting solutions visualised using the functions tree_plot (provides a visualisation of the entire model) and node_plot (provides a more detailed visualisation of a single node in the hierarchical model). Nodes can then be removed using the function tree_prune, or added using the function tree_split, depending on the apparent validity of the existing clustering model.
}
\details{

Package:  PPCI

Type: Package

Title:  Projection Pursuit for Cluster Identification

Version:  0.1.0

Depends:  R (>= 2.10.0)

License:  GPL-3

LazyData: yes

Imports:  rARPACK

}
\author{
David Hofmeyr[aut, cre] and Nicos Pavlidis[aut]

Maintainer: \packageMaintainer{PPCI}
}
\seealso{A \bold{MATLAB} toolbox with the same functionality as this package is available at \url{https://github.com/nicospavlidis/PPCI-MATLAB}. Outputs may differ slightly due to differences between \bold{R} and \bold{MATLAB}'s base optimisation software.}
\references{
Pavlidis N.G., Hofmeyr D.P., Tasoulis S.K. (2016) Minimum Density Hyperplanes. \emph{Journal of
Machine Learning Research}, \bold{17}(156), 1--33.

Hofmeyr, D., Pavlidis, N. (2015) Maximum Clusterability Divisive Clustering. \emph{Computational Intelligence, 2015 IEEE Symposium Series on}, pp. 780--786.

Hofmeyr, D. (2016) Clustering by Minimum Cut Hyperplanes. \emph{IEEE Transactions on Pattern Analysis and Machine Intelligence}, \bold{39}(8), 1547 -- 1560.
}
\keyword{ package }
