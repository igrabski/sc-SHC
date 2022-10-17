# Significance Analysis for Clustering Single-Cell RNA-Sequencing Data

We introduce a model-based hypothesis testing approach for evaluating single-cell RNA-sequencing (scRNA-seq) clusters. This approach is implemented in two ways: (1) a stand-alone clustering pipeline with built-in hypothesis testing to produce clusters corresponding to distinct cell populations and (2) a post-hoc method that can evaluate the statistical significance of any provided set of clusters. 

# Usage

The functions in `significance_analysis.R` represent those used in our preprint. We also now offer a faster version in `significance_analysis_fast.R`, which makes some approximations to the method; please note that this is a developmental version. 

To use the stand-alone clustering pipeline, sc-SHC (single-cell significance of hierarchical clustering), the following command can be used:

```
# For the original version
source('significance_analysis.R')

clusters <- scSHC(data)

# For the developmental, faster version
source('significance_analysis_fast.R')

clusters <- scSHC(data)
```

Here, ```data``` should be a matrix where the rows are genes and the columns are cells. Optionally, the following parameters can be adjusted: 

* ```alpha```, which controls the family-wise error rate (default 0.05). If the goal is discovery, consider setting a more lenient alpha, such as 0.25.
* ```num_PCs```, which controls the number of principal components (default 30).
* ```num_features```, which controls the number of genes used (default 2500). 

To evaluate the significance of any provided set of clusters, the following command can be used:

```
# For the original version
source('significance_analysis.R')

new_clusters <- testClusters(data, as.character(clusters))

# For the developmental, faster version
source('significance_analysis_fast.R')

new_clusters <- testClusters(data, as.character(clusters))
```

Here, ```data``` is the same as before, and ```clusters``` should be a vector of cluster labels, corresponding to cells in the same order as the columns of the ```data``` matrix. The same parameters as above can be adjusted. Additionally, if desired, a given set of genes can be provided through the parameter ```var.genes``` rather than allowing our approach to identify informative genes on its own.
