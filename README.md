# Significance Analysis for Clustering Single-Cell RNA-Sequencing Data

Contact information: igrabski[at]nygenome[dot]org

We introduce a model-based hypothesis testing approach for evaluating single-cell RNA-sequencing (scRNA-seq) clusters. This approach is implemented in two ways: (1) a stand-alone clustering pipeline with built-in hypothesis testing to produce clusters corresponding to distinct cell populations and (2) a post-hoc method that can evaluate the statistical significance of any provided set of clusters. 

Our package can be installed as follows:

```
# install.packages("devtools")
devtools::install_github("igrabski/sc-SHC")
```

# Usage

To use the stand-alone clustering pipeline, sc-SHC (single-cell significance of hierarchical clustering), the following command can be used:

```
library(scSHC)
clusters <- scSHC(data)
```

Here, ```data``` should be a (possibly sparse) matrix where the rows are genes and the columns are cells. Optionally, the following parameters can be adjusted: 

* ```batch```, which is `NULL` by default, can be a character vector of known batch labels.
* ```alpha```, which controls the family-wise error rate (default 0.05). If the goal is discovery, consider setting a more lenient alpha, such as 0.25.
* ```num_features```, which controls the number of genes used (default 2500). 
* ```num_PCs```, which controls the number of principal components (default 30).
* ```parallel```, which is `TRUE` by default, can be set to `FALSE` to disable parallelization.
* ```cores```, which controls the number of cores used if `parallel = T` (default 2).

To evaluate the significance of any provided set of clusters, the following command can be used:

```
library(scSHC)
new_clusters <- testClusters(data, as.character(clusters))
```

Here, ```data``` is the same as before, and ```clusters``` should be a character vector of cluster labels, corresponding to cells in the same order as the columns of the ```data``` matrix. The same parameters as above can be adjusted. Additionally, if desired, a given set of genes can be provided through the parameter ```var.genes``` rather than allowing our approach to identify informative genes on its own.
