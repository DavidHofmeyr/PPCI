\name{plot.ppci_projection_solution}
\alias{plot.ppci_projection_solution}
\title{Visualise a Data Set Projected from Projection Pursuit}
\description{
  Provides a visualisation of a data set using a low dimensional projection obtained using projection pursuit based on cluster based projection indices.
}
\usage{
\method{plot}{ppci_projection_solution}(x, labels = NULL, pairs = NULL, PCA = NULL, transparency = NULL, ...)
}
\arguments{
  \item{x}{a solution arising from one of the functions mcdr, mddr and ncutdr.}
  \item{labels}{(optional) a vector of class labels. if provided then points in different classes are plotted in different colours.}
  \item{pairs}{(optional) if omitted then the first two dimensions are shown. If an integer > 1 then all pairs from the first 1:pairs dimensions are shown.}
  \item{PCA}{(optional) logical. If TRUE then an additional principal component rotation is applied. This can make more clusters visible within a low dimensional projection.}
  \item{transparency}{(optional) if ommitted then points in scatterplots are
  shown as solid. If set to a value in (0, 1) then points are shown with transparency. Lower values correspond with a greater degree of transparency.}
  \item{...}{additional graphical parameters. Currently none implemented.}
}
\seealso{mcdr, mddr, ncutdr}
\value{None}
\keyword{internal}
