\name{testClusters}
\alias{testClusters}
\title{Significance Analysis for Pre-Computed Clusters}
\usage{
testClusters(data, cluster_ids, batch = NULL, var.genes = NULL, alpha = 0.05,
    num_features = 2500, num_PCs = 30, parallel = T, cores = 2)
}
\description{
Performs significance analysis on pre-computed clusters in single-cell data.
}
\arguments{
\item{data}{raw counts \code{matrix} or \code{Matrix} in genes by cells format,
with row names and column names}

\item{cluster_ids}{character vector of pre-computed cluster ids, in the same
order as the columns of \code{data}}

\item{batch}{character vector of batch labels, in the same order as the
columns of \code{data}}

\item{var.genes}{(optional) vector of gene names representing the subset of
genes used for clustering}

\item{alpha}{family-wise error rate}

\item{num_features}{number of top genes to include in analysis; not used if
\code{var.genes} is supplied}

\item{num_PCs}{number of PCs to use in analysis}

\item{parallel}{whether or not parallelization should be used}

\item{cores}{number of cores; ignored if \code{parallel = F}}
}
\value{
  A list containing the vector of final cluster labels, and a \code{data.tree}
  object of the adjusted p-values for every split considered in the hierarchical
  clustering tree.
}
\examples{
data(counts)
cluster_labels <- c(rep('cluster1',288),rep('cluster2',289),
                    rep('cluster3',288),rep('cluster4',289))
final_labels <- testClusters(counts, cluster_labels)
}
