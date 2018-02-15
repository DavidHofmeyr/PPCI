\name{cluster_performance}
\alias{cluster_performance}
\title{External Cluster Validity Metrics}
\description{
  Computes four popular external cluster validity metrics (adjusted Rand index, purity, V-measure and Normalised Mutual Information) through comparison of cluster assignments and true class labels. 
}
\usage{
cluster_performance(assigned, labels, beta)
}
\arguments{
  \item{assigned}{a vector of cluster assignments made by a clustering algorithm.}
  \item{labels}{a vector of true class labels to be compared with assigned.}
  \item{beta}{(optional) positive numeric, used in the computation of V-measure. larger values apply higher weight to homogeneity over completeness measures. if omitted then beta = 1 (equal weight applied to both measures).}
}
\value{
  a vector containing the four evaluation metrics listed in the description.
}
\references{
Zhao Y., Karypis G. (2004) Empirical and Theoretical Comparisons of Selected Criterion
Functions for Document Clustering. \emph{Machine Learning}, \bold{55}(3), 311--331.

Strehl A., Ghosh J. (2002) Cluster ensembles—a knowledge reuse framework for combining
multiple partitions. \emph{Journal of Machine Learning Research}, \bold{3}, 583--617.

Rosenberg A., Hirschberg J. (2007) V-Measure: A Conditional Entropy-Based External
Cluster Evaluation Measure. \emph{EMNLP-CoNLL}, \bold{7}, 410--420. Citeseer.

Hubert, L., Arabie, P. (1985) Comparing Partitions. \emph{Journal of Classification}, \bold{2}(1), 193--218.
}
\examples{
## load dermatology dataset
data(dermatology)

## obtain clustering solution using MCDC
sol <- mcdc(dermatology$x, 6)

## evaluate solution using external cluster validity measures
cluster_performance(sol$cluster, dermatology$c)
}
\keyword{file}
