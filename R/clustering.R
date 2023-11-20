#' Fit model for one batch
#'
#' @importFrom Matrix t colMeans rowMeans
#' @importFrom RSpectra eigs_sym
#' @importFrom sfsmisc posdefify
#'
#' @noRd
fit_model_batch <- function(y,on_genes,num_PCs) {
  # Compute sample moments of the on genes
  on_counts <- Matrix::t(y[on_genes,])
  cov <- cov(as.matrix(on_counts))
  means <- Matrix::colMeans(on_counts)

  # Use method of moments to estimate parameters
  sigmas <- log(((diag(cov)-means)/means^2)+1)
  mus <- log(means)-0.5*sigmas
  mus.sum <- tcrossprod(array(mus,dim=c(length(mus),1)),
                        array(1,dim=c(length(mus),1)))+
    tcrossprod(array(1,dim=c(length(mus),1)),array(mus,dim=c(length(mus),1)))
  sigmas.sum <- tcrossprod(array(sigmas,dim=c(length(sigmas),1)),
                           array(1,dim=c(length(sigmas),1)))+
    tcrossprod(array(1,dim=c(length(sigmas),1)),
               array(sigmas,dim=c(length(sigmas),1)))
  rhos <- suppressWarnings(log(cov/(exp(mus.sum+0.5*sigmas.sum))+1))
  rhos[is.na(rhos)|is.infinite(rhos)] <- (-10)
  diag(rhos) <- sigmas

  # Make the covariance matrix positive-definite
  on_cov_eigs <- RSpectra::eigs_sym(as.matrix(rhos),k=min(c(nrow(rhos)-1,num_PCs)))
  num_pos <- sum(on_cov_eigs$values>0)
  on_cov.sub <- on_cov_eigs$vectors[,1:num_pos]%*%
    sqrt(diag(on_cov_eigs$values[1:num_pos]))
  on_cov <- tcrossprod(on_cov.sub)
  diag(on_cov) <- diag(rhos)
  on_cov <- sfsmisc::posdefify(on_cov)
  on_cov.sqrt <- t(chol(on_cov))

  return(list(Matrix::rowMeans(y),mus,on_cov.sqrt))
}

#' Fit model
#'
#' @noRd
#'
fit_model <- function(y,on_genes,x,num_PCs) {
  on_means <- list()
  on_cov.sqrt <- list()
  lambdas <- list()

  for (b in unique(x)) {
    params <- fit_model_batch(y[,x==b],on_genes=on_genes,num_PCs=num_PCs)
    lambdas[[as.character(b)]] <- params[[1]]
    on_means[[as.character(b)]] <- params[[2]]
    on_cov.sqrt[[as.character(b)]] <- params[[3]]
  }

  return(list(lambdas,on_means,on_cov.sqrt))
}

#' Generate one null sample
#'
#' @importFrom stats rpois rnorm
#'
#' @noRd
#'
generate_null <- function(y,params,on_genes,x) {
  lambdas <- params[[1]]
  on_means <- params[[2]]
  on_cov.sqrt <- params[[3]]

  null <- array(0,dim=dim(y))
  rownames(null) <- rownames(y)

  for (b in unique(x)) {
    num_gen <- min(sum(x==b),1000)
    names(lambdas[[as.character(b)]]) <- rownames(y)
    null[-on_genes,which(x==b)[1:num_gen]] <-
      array(rpois(num_gen*(nrow(null)-length(on_genes)),
                  lambdas[[as.character(b)]][-on_genes]),
            dim=c(nrow(null)-length(on_genes),num_gen))
    Y <- exp(sweep(on_cov.sqrt[[as.character(b)]]%*%
                     array(rnorm(num_gen*length(on_genes)),
                           dim=c(length(on_genes),num_gen)),1,
                   on_means[[as.character(b)]],'+'))
    null[on_genes,which(x==b)[1:num_gen]] <-
      array(rpois(length(Y),Y),dim=dim(Y))
  }

  return(list(null[,colSums(null)>0],x[colSums(null)>0]))
}

#' Generate one null sample and compute test statistic
#'
#' @importFrom fastcluster hclust
#' @importFrom stats dist
#' @importFrom BiocNeighbors queryKNN AnnoyParam
#'
#' @noRd
#'
generate_null_statistic <- function(y,params,on_genes,x,num_PCs,
                                    gm,labs,posthoc) {
  null_set <- generate_null(y,params,on_genes,x)
  null <- null_set[[1]]
  batch_null <- null_set[[2]]

  if (!posthoc) {
    null_gm <- reduce_dimension(null,batch_null,num_PCs)[[2]]
    null_gm.d <- dist(null_gm)
    hc2 <- cutree(hclust(null_gm.d,method='ward.D'),2)
  } else {
    null_gm <- reduce_dimension(null,batch_null,num_PCs)[[2]]
    pdev <- poisson_dev_batch(null,batch_null)
    pdev <- t(scale(Matrix::t(pdev),scale=F))
    gm2 <- t(crossprod(gm[[1]]$vectors,pdev))
    knns <- queryKNN(gm[[2]], gm2,
                     k=15, BNPARAM=AnnoyParam())$index
    neighbor.labels <- apply(knns,1,function(x) labs[x])
    hc2 <- unlist(apply(neighbor.labels,2,function(x)
      sample(names(table(x))[as.vector(table(x))==max(table(x))],1)))

    if (length(unique(hc2))==1) {
      null_gm.d <- dist(null_gm)
      hc2 <- cutree(hclust(null_gm.d,method='ward.D'),2)
    }
  }

  Qclust2 <- sapply(unique(batch_null),function(b)
    if (length(unique(hc2[batch_null==b]))==2 &
        min(table(hc2[batch_null==b]))>=2) {
      ward_linkage(null_gm[batch_null==b,],hc2[batch_null==b])
    } else {
      0
    })
  names(Qclust2) <- unique(batch_null)
  return(median(Qclust2))
}

#' Test one split
#'
#' @importFrom stats pnorm median
#' @importFrom parallel mclapply
#' @importFrom MASS fitdistr
#' @importFrom matrixStats rowMins
#'
#' @noRd
#'
test_split <- function(data,ids1,ids2,var.genes,num_PCs,batch,
                       alpha_level,cores,posthoc) {
  # Re-order data
  cell1s <- data[,ids1]
  cell2s <- data[,ids2]
  true <- cbind(cell1s,cell2s)
  batch <- batch[c(ids1,ids2)]

  # Re-run dimension reduction and calculate test statistic
  gm <- reduce_dimension(true[var.genes,],batch,num_PCs)
  gm_sub.x <- gm[[2]]
  labs <- c(rep(1,length(ids1)),rep(2,length(ids2)))
  Qclust <- sapply(unique(batch),function(b)
    ward_linkage(gm_sub.x[batch==b,],labs[batch==b]))
  names(Qclust) <- unique(batch)
  stat <- median(Qclust)

  # Determine the set of "on" genes
  phi_stat <- poisson_dispersion_stats(true[var.genes,])
  check_means <- matrixStats::rowMins(sapply(unique(batch),function(b)
    Matrix::rowSums(true[var.genes,batch==b])))
  on_genes <- which(pnorm(phi_stat,lower.tail=F)<0.05&check_means!=0)

  # Fit model
  params <- fit_model(true[var.genes,],on_genes,batch,num_PCs)

  # Generate null distribution of test statistics
  Qclusts2_1 <- mclapply(1:10,function(i) {
    generate_null_statistic(true[var.genes,],params,on_genes,batch,
                            num_PCs,gm,labs,posthoc)
  },mc.cores=cores)

  # Quit early if p-value is much smaller or much larger than alpha level
  Qclusts2 <- unlist(Qclusts2_1)
  fit <- fitdistr(Qclusts2,'normal')
  pval <- 1-pnorm(stat,mean=fit$estimate[1],sd=fit$estimate[2])
  if (pval < 0.1*alpha_level | pval > 10*alpha_level) {
    return(pval)
  }

  # Otherwise, keep going
  Qclusts2_2 <- mclapply(11:50,function(i) {
    generate_null_statistic(true[var.genes,],params,on_genes,
                            batch,num_PCs,gm,labs,posthoc)
  },mc.cores=cores)

  # Compute smoothed p-value
  Qclusts2 <- c(Qclusts2,unlist(Qclusts2_2))
  fit <- fitdistr(Qclusts2,'normal')
  return(1-pnorm(stat,mean=fit$estimate[1],sd=fit$estimate[2]))
}

#' scSHC
#'
#' Significance of Hierarchical Clustering for Single-Cell Data
#'
#' @details
#' Full clustering pipeline with built-in hypothesis testing.
#' Performs hierarchical clustering on single-cell data, with significance
#' analysis built into the algorithm.
#'
#' @param data raw counts `matrix` or `Matrix` in genes by cells format,
#' with row names and column names
#' @param batch character vector of batch labels, in the same order as the
#' columns of `data`
#' @param alpha family-wise error rate
#' @param num_features number of top genes to include in analysis
#' @param num_PCs number of PCs to use in analysis
#' @param parallel whether or not parallelization should be used
#' @param cores number of cores; ignored if `parallel = FALSE`
#'
#' @return A list containing the vector of final cluster labels, and a
#' `data.tree` object of the adjusted p-values for every split considered
#' in the hierarchical clustering tree.
#'
#' @export
#'
#' @importFrom dendextend cutree get_leaves_attr
#' @importFrom scry devianceFeatureSelection
#' @importFrom data.tree Node
#' @importFrom matrixStats rowMins
#'
#' @examples
#' data(counts)
#' \dontrun{
#'   clusters <- scSHC(counts)
#' }
scSHC <- function(data,batch=NULL,alpha=0.05,num_features=2500,
                  num_PCs=30,parallel=T,cores=2) {
  if (!parallel) {
    cores <- 1
  }
  if (is.null(batch)) {
    batch <- rep("1",ncol(data))
  }
  if (is.factor(batch)|is.numeric(batch)) {
    warning("Converting batch vector to character")
    batch <- as.character(batch)
  }
  names(batch) <- colnames(data)
  if (is.null(colnames(data))) {
    warning("Assigning cell names automatically to the columns")
    colnames(data) <- paste0('cell',1:ncol(data))
  }

  # Get variable features
  dev <- scry::devianceFeatureSelection(data)
  var.genes <- rownames(data)[order(dev,decreasing=T)[1:num_features]]

  # Dimension reduction and clustering
  gm.x <- reduce_dimension(data[var.genes,],batch,num_PCs)[[2]]
  gm.d <- dist(gm.x)
  hcl <- fastcluster::hclust(gm.d,method='ward.D')
  dend <- as.dendrogram(hcl)

  # Traverse the tree and store clusters and q-FWERs
  dends_to_test <- list(dend)
  clusters <- list()
  node0 <- NULL
  counter <- 0
  parents <- list('root')

  while (length(dends_to_test)>0) {
    # Identify the two-way split
    cuts <- dendextend::cutree(dends_to_test[[1]],k=2,order_clusters_as_data=F)
    leaves <- get_leaves_attr(dends_to_test[[1]],'label')
    ids1 <- leaves[cuts==1]
    ids2 <- leaves[cuts==2]
    alpha_level <- alpha*((length(leaves)-1)/(ncol(data)-1))

    # Remove any cells from batches with poor representation
    tab <- table(batch[c(ids1,ids2)],cuts)
    to.keep <- rownames(tab)[which(matrixStats::rowMins(tab)>20)]
    ids1 <- ids1[batch[ids1]%in%to.keep]
    ids2 <- ids2[batch[ids2]%in%to.keep]

    # Get p-value of the split
    if (length(to.keep)>0) {
      test <- test_split(data,ids1,ids2,var.genes,num_PCs,
                         batch,alpha_level,cores,posthoc=F)
    } else {
      test <- 1
    }

    # Compare to significance threshold
    if (test < alpha_level) {
      # If significant: add left and right branches to testing stack
      left <- dends_to_test[[1]][[1]]
      right <- dends_to_test[[1]][[2]]
      dends_to_test[[length(dends_to_test)+1]] <- left
      dends_to_test[[length(dends_to_test)+1]] <- right

      # Compute q-FWER
      if (is.null(node0)) {
        node0 <- Node$new(paste0('Node 0: ',
           min(round(test*(ncol(data)-1)/(length(leaves)-1),2),1)))
      } else {
        do.call("<-",list(paste0('node',counter),
           eval(parse(text=parents[[1]]))$AddChild(paste0('Node ',counter,': ',
           min(round(test*(ncol(data)-1)/(length(leaves)-1),2),1)))))
      }
      parents[[length(parents)+1]] <- paste0('node',counter)
      parents[[length(parents)+1]] <- paste0('node',counter)
      counter <- counter + 1
    } else {
      # If not significant: create cluster
      clusters[[length(clusters)+1]] <- leaves

      # Compute q-FWER
      if (is.null(node0)) {
        node0 <- Node$new(paste0('Cluster 0: ',
          min(round(test*(ncol(data)-1)/(length(leaves)-1),2),1)))
      } else {
        do.call("<-",list(paste0('Cluster',length(clusters)),
          eval(parse(text=parents[[1]]))$AddChild(paste0('Cluster ',length(clusters),': ',
          min(round(test*(ncol(data)-1)/(length(leaves)-1),2),1)))))
      }
    }
    dends_to_test[[1]] <- NULL
    parents[[1]] <- NULL
  }

  # Produce vector of cluster labels
  cluster_labels <- rep(0,ncol(data))
  names(cluster_labels) <- colnames(data)
  for (i in 1:length(clusters)) {
    cluster_labels[clusters[[i]]] <- i
  }

  return(list(cluster_labels,node0))
}

#' testClusters
#'
#' Significance Analysis for Pre-Computed Clusters
#'
#' @details
#' Performs significance analysis on pre-computed clusters in single-cell data.
#'
#' @param data raw counts `matrix` or `Matrix` in genes by cells format,
#' with row names and column names
#' @param cluster_ids character vector of pre-computed cluster ids, in the same
#' order as the columns of `data`
#' @param batch character vector of batch labels, in the same order as the
#' columns of `data`
#' @param var.genes (optional) vector of gene names representing the subset of
#' genes used for clustering
#' @param alpha family-wise error rate
#' @param num_features number of top genes to include in analysis; not used if
#' `var.genes` is supplied
#' @param num_PCs number of PCs to use in analysis
#' @param parallel whether or not parallelization should be used
#' @param cores number of cores; ignored if `parallel = FALSE`
#'
#' @return A list containing the vector of final cluster labels, and a
#' `data.tree` object of the adjusted p-values for every split considered in
#' the hierarchical clustering tree.
#'
#' @export
#'
#' @importFrom stats aggregate as.dendrogram
#' @importFrom scry devianceFeatureSelection
#' @importFrom dendextend get_branches_heights
#' @importFrom matrixStats rowMins
#'
#' @examples
#' data(counts)
#' cluster_labels <- c(rep('cluster1',288),rep('cluster2',289),
#'                     rep('cluster3',288),rep('cluster4',289))
#' \dontrun{
#'   final_labels <- testClusters(counts, cluster_labels)
#' }
testClusters <- function(data,cluster_ids,batch=NULL,var.genes=NULL,
                         alpha=0.05,num_features=2500,num_PCs=30,
                         parallel=T,cores=2) {
  if (is.null(var.genes)) {
    dev <- scry::devianceFeatureSelection(data)
    var.genes <- rownames(data)[order(dev,decreasing=T)[1:num_features]]
  } else {
    var.genes <- intersect(var.genes,rownames(data))
  }
  if (is.null(batch)) {
    batch <- rep("1",ncol(data))
  }
  if (is.factor(batch)|is.numeric(batch)) {
    warning("Converting batch vector to character")
    batch <- as.character(batch)
  }
  names(batch) <- colnames(data)
  if (is.factor(cluster_ids)|is.numeric(cluster_ids)) {
    warning("Converting cluster labels to character")
    cluster_ids <- as.character(cluster_ids)
  }
  if (is.null(colnames(data))) {
    warning("Assigning cell names automatically to the columns")
    colnames(data) <- paste0('cell',1:ncol(data))
  }

  # Create hierarchy of clusters
  pseudobulk <- aggregate(t(as.matrix(data[var.genes,])),
                          by=list(cluster_ids),sum)
  rownames(pseudobulk) <- pseudobulk[,1]
  pseudobulk <- pseudobulk[,2:ncol(pseudobulk)]
  pseudobulk <- sweep(pseudobulk,1,rowSums(pseudobulk),'/')
  pseudobulk.d <- dist(pseudobulk)
  pseudobulk.hc <- hclust(pseudobulk.d,method='ward.D')
  dend <- as.dendrogram(pseudobulk.hc)

  # Traverse the tree and store clusters and q-FWERs
  dends_to_test <- list(dend)
  clusters <- list()
  node0 <- NULL
  counter <- 0
  parents <- list('root')

  while (length(dends_to_test)>0) {
    # Identify the two-way split
    cuts <- dendextend::cutree(dends_to_test[[1]],k=2,order_clusters_as_data=F)
    leaves <- get_leaves_attr(dends_to_test[[1]],'label')
    ids1 <- which(cluster_ids%in%leaves[cuts==1])
    ids2 <- which(cluster_ids%in%leaves[cuts==2])
    alpha_level <- alpha*((sum(cluster_ids%in%leaves)-1)/(ncol(data)-1))

    # Remove any cells from batches with poor representation
    tab <- table(batch[c(ids1,ids2)],c(rep(1,length(ids1)),rep(2,length(ids2))))
    to.keep <- rownames(tab)[which(matrixStats::rowMins(tab)>20)]
    ids1 <- ids1[batch[ids1]%in%to.keep]
    ids2 <- ids2[batch[ids2]%in%to.keep]

    # Get p-value of the split
    if (length(to.keep)>0) {
      test <- test_split(data,ids1,ids2,var.genes,num_PCs,batch,
                         alpha_level,cores,posthoc=T)
    } else {
      test <- 1
    }

    # Compare to significance threshold
    if (test < alpha_level) {
      # If significant: get left and right branches
      left <- dends_to_test[[1]][[1]]
      right <- dends_to_test[[1]][[2]]

      # Compute q-FWER
      if (is.null(node0)) {
        node0 <- Node$new(paste0('Node 0: ',
          min(round(test*(ncol(data)-1)/(sum(cluster_ids%in%leaves)-1),2),1)))
      } else {
        do.call("<-",list(paste0('node',counter),
         eval(parse(text=parents[[1]]))$AddChild(paste0('Node ',counter,': ',
         min(round(test*(ncol(data)-1)/(sum(cluster_ids%in%leaves)-1),2),1)))))
      }

      # If there are no further splits down a branch, make it a cluster
      # Otherwise, add it to the stack
      if (suppressWarnings(max(get_branches_heights(left)))!=(-Inf)) {
        dends_to_test[[length(dends_to_test)+1]] <- left
        parents[[length(parents)+1]] <- paste0('node',counter)
      } else {
        clusters[[length(clusters)+1]] <- get_leaves_attr(left,'label')
        do.call("<-",list(paste0('Cluster',length(clusters)),
          eval(parse(text=paste0('node',counter)))$AddChild(paste0('Cluster ',
          length(clusters)))))
      }
      if (suppressWarnings(max(get_branches_heights(right)))!=(-Inf)) {
        dends_to_test[[length(dends_to_test)+1]] <- right
        parents[[length(parents)+1]] <- paste0('node',counter)
      } else {
        clusters[[length(clusters)+1]] <- get_leaves_attr(right,'label')
        do.call("<-",list(paste0('Cluster',length(clusters)),
          eval(parse(text=paste0('node',counter)))$AddChild(paste0('Cluster ',
          length(clusters)))))
      }

      counter <- counter + 1
    } else {
      # If not significant: create cluster
      clusters[[length(clusters)+1]] <- leaves

      # Compute q-FWER
      if (is.null(node0)) {
        node0 <- Node$new(paste0('Cluster 0: ',
          min(round(test*(ncol(data)-1)/(sum(cluster_ids%in%leaves)-1),2),1)))
      } else {
        do.call("<-",list(paste0('Cluster',length(clusters)),
         eval(parse(text=parents[[1]]))$AddChild(paste0('Cluster ',length(clusters),': ',
         min(round(test*(ncol(data)-1)/(sum(cluster_ids%in%leaves)-1),2),1)))))
      }
    }
    dends_to_test[[1]] <- NULL
    parents[[1]] <- NULL
  }

  # Return vector of new, merged cluster labels
  cluster_labels <- cluster_ids
  for (x in 1:length(clusters)) {
    cluster_labels[which(cluster_ids%in%clusters[[x]])] <- paste0('new',x)
  }
  return(list(cluster_labels,node0))
}


#' A sample count matrix
#'
#' @name counts
#' @docType data
NULL

#' A(nother) sample count matrix, just stored in bzip2 compressed format
#'
#' @name counts_compressed
#' @docType data
NULL
