\name{scSHC}
\alias{scSHC}
\title{Significance of Hierarchical Clustering for Single-Cell Data}
\usage{
scSHC(data, batch = NULL, alpha = 0.05, num_features = 2500, num_PCs = 30,
    parallel = T, cores = 2)
}
\description{
Performs hierarchical clustering on single-cell data, with significance analysis
built into the algorithm.
}
\arguments{
\item{data}{raw counts \code{matrix} or \code{Matrix} in genes by cells format,
with row names and column names}

\item{batch}{character vector of batch labels, in the same order as the
columns of \code{data}}

\item{alpha}{family-wise error rate}

\item{num_features}{number of top genes to include in analysis}

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
clusters <- scSHC(counts)
}
