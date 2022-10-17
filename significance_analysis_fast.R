library(scry)
library(AER)
library(MASS)
library(cluster)
library(sctransform)
library(irlba)
library(dendextend)
library(BiocNeighbors)
library(parallel)

## Test significance of a split in a hierarchical tree as part of sc-SHC (parallelized).
# data = Genes x cells UMI matrix.
# ids1 = Cell indices belonging to the left branch.
# ids2 = Cell indices belonging to the right branch.
# var.genes = Informative gene names.
# num_PCs = Number of principal components to use.
# cores = Number of cores to use.

test_split.parallel <- function(data,ids1,ids2,var.genes,num_PCs,cores,alpha_level) { 
  # Re-order data
  cell1s <- data[,ids1]
  cell2s <- data[,ids2]
  true <- cbind(cell1s,cell2s)
  
  # Re-run dimension reduction on just these clusters
  pdev <- poisson.dev(true[var.genes,])
  pdev <- t(scale(t(pdev),scale=F))
  gm_sub <- RSpectra::eigs_sym(as.matrix(tcrossprod(pdev)),k=num_PCs)
  gm_sub.x <- t(crossprod(gm_sub$vectors,pdev))
  gm_sub.d <- dist(gm_sub.x[,1:2])
  
  # Compute average silhouette
  if (length(ids1)>1&length(ids2)>1) {
    labs <- c(rep(1,length(ids1)),rep(2,length(ids2)))
    s <- silhouette(as.integer(labs),dist=gm_sub.d)
    Qclust <- mean(s[,3])
  } else {
    Qclust <- mean(gm_sub.d)
  }
  
  # Determine the set of "on" genes
  phi_stat <- poisson_dispersion_stats(true[var.genes,])
  on_genes <- which(pnorm(phi_stat,lower.tail=F)<0.05)
  
  # Compute sample moments of the on genes
  on_counts <- t(true[var.genes[on_genes],])
  cov <- cov(as.matrix(on_counts))
  means <- colMeans(on_counts)
  
  # Use method of moments to estimate parameters
  sigmas <- log(((diag(cov)-means)/means^2)+1)
  mus <- log(means)-0.5*sigmas
  mus.sum <- tcrossprod(array(mus,dim=c(length(mus),1)),array(1,dim=c(length(mus),1)))+
    tcrossprod(array(1,dim=c(length(mus),1)),array(mus,dim=c(length(mus),1)))
  sigmas.sum <- tcrossprod(array(sigmas,dim=c(length(sigmas),1)),array(1,dim=c(length(sigmas),1)))+
    tcrossprod(array(1,dim=c(length(sigmas),1)),array(sigmas,dim=c(length(sigmas),1)))
  rhos <- log(cov/(exp(mus.sum+0.5*sigmas.sum))+1)
  rhos[is.na(rhos)|is.infinite(rhos)] <- (-10)
  diag(rhos) <- sigmas
  
  on_cov_eigs <- RSpectra::eigs_sym(as.matrix(rhos),k=num_PCs)
  num_pos <- sum(on_cov_eigs$values>0)
  on_cov.sub <- on_cov_eigs$vectors[,1:num_pos]%*%sqrt(diag(on_cov_eigs$values[1:num_pos]))
  on_cov <- tcrossprod(on_cov.sub)
  diag(on_cov) <- diag(rhos)
  on_cov <- sfsmisc::posdefify(on_cov)
  on_cov.sqrt <- t(chol(on_cov))
  on_means <- mus
  lambdas <- rowMeans(true[var.genes,])
  
  nm <- min(5000,ncol(true))
  
  # Generate null distribution of test statistics
  if (nm==5000) {
    Qclusts2_1 <- mclapply(1:10,function(i) {
      # Generate a null dataset under the model
      null <- array(0,dim=c(length(var.genes),nm))
      rownames(null) <- names(lambdas) <- var.genes
      null[-on_genes,] <- array(rpois(ncol(null)*length(var.genes[-on_genes]),lambdas[-on_genes]),
                                dim=c(length(var.genes[-on_genes]),ncol(null)))
      Y <- exp(sweep(on_cov.sqrt%*%array(rnorm(nm*length(on_genes)),
                                         dim=c(length(on_genes),nm)),1,on_means,'+'))
      null[on_genes,] <- array(rpois(length(Y),Y),dim=dim(Y))
      
      # Reduce dimensionality of null dataset and cluster
      pdev2 <- poisson.dev(null)
      pdev2 <- t(scale(t(pdev2),scale=F))
      gm2 <- RSpectra::eigs_sym(as.matrix(tcrossprod(pdev2)),k=num_PCs)
      gm2.x <- t(crossprod(gm2$vectors,pdev2))
      gm2.d <- dist(gm2.x)
      hc2 <- cutree(hclust(gm2.d,method='ward.D2'),2)
      
      # Compute average silhouette
      if (length(ids1)>1&length(ids2)>2) {
        s <- silhouette(as.integer(hc2),dist=dist(gm2.x[,1:2]))
        mean(s[,3])
      } else {
        mean(gm2.d)
      }
    },mc.cores=cores)
    
    # Quit early if p-value is much smaller than alpha level
    Qclusts2 <- unlist(Qclusts2_1)
    fit <- fitdistr(Qclusts2,'normal')
    pval <- 1-pnorm(Qclust,mean=fit$estimate[1],sd=fit$estimate[2])
    if (pval < 0.1*alpha_level) {
      return(pval)
    }
    
    # Otherwise, keep going
    Qclusts2_2 <- mclapply(11:50,function(i) {
      # Generate a null dataset under the model
      null <- array(0,dim=c(length(var.genes),nm))
      rownames(null) <- names(lambdas) <- var.genes
      null[-on_genes,] <- array(rpois(ncol(null)*length(var.genes[-on_genes]),lambdas[-on_genes]),
                                dim=c(length(var.genes[-on_genes]),ncol(null)))
      Y <- exp(sweep(on_cov.sqrt%*%array(rnorm(nm*length(on_genes)),
                                         dim=c(length(on_genes),nm)),1,on_means,'+'))
      null[on_genes,] <- array(rpois(length(Y),Y),dim=dim(Y))
      
      # Reduce dimensionality of null dataset and cluster
      pdev2 <- poisson.dev(null)
      pdev2 <- t(scale(t(pdev2),scale=F))
      gm2 <- RSpectra::eigs_sym(as.matrix(tcrossprod(pdev2)),k=num_PCs)
      gm2.x <- t(crossprod(gm2$vectors,pdev2))
      gm2.d <- dist(gm2.x)
      hc2 <- cutree(hclust(gm2.d,method='ward.D2'),2)
      
      # Compute average silhouette
      if (length(ids1)>1&length(ids2)>2) {
        s <- silhouette(as.integer(hc2),dist=dist(gm2.x[,1:2]))
        mean(s[,3])
      } else {
        mean(gm2.d)
      }
    },mc.cores=cores)
    
    # Compute smoothed p-value
    Qclusts2 <- c(Qclusts2,unlist(Qclusts2_2))
    fit <- fitdistr(Qclusts2,'normal')
    return(1-pnorm(Qclust,mean=fit$estimate[1],sd=fit$estimate[2]))
  } else {
    Qclusts2 <- mclapply(1:50,function(i) {
      # Generate a null dataset under the model
      null <- array(0,dim=c(length(var.genes),nm))
      rownames(null) <- names(lambdas) <- var.genes
      null[-on_genes,] <- array(rpois(ncol(null)*length(var.genes[-on_genes]),lambdas[-on_genes]),
                                dim=c(length(var.genes[-on_genes]),ncol(null)))
      Y <- exp(sweep(on_cov.sqrt%*%array(rnorm(nm*length(on_genes)),
                                         dim=c(length(on_genes),nm)),1,on_means,'+'))
      null[on_genes,] <- array(rpois(length(Y),Y),dim=dim(Y))
      
      # Reduce dimensionality of null dataset and cluster
      pdev2 <- poisson.dev(null)
      pdev2 <- t(scale(t(pdev2),scale=F))
      gm2 <- RSpectra::eigs_sym(as.matrix(tcrossprod(pdev2)),k=num_PCs)
      gm2.x <- t(crossprod(gm2$vectors,pdev2))
      gm2.d <- dist(gm2.x)
      hc2 <- cutree(hclust(gm2.d,method='ward.D2'),2)
      
      # Compute average silhouette
      if (length(ids1)>1&length(ids2)>2) {
        s <- silhouette(as.integer(hc2),dist=dist(gm2.x[,1:2]))
        mean(s[,3])
      } else {
        mean(gm2.d)
      }
    },mc.cores=cores)
    
    # Compute smoothed p-value
    Qclusts2 <- unlist(Qclusts2)
    fit <- fitdistr(Qclusts2,'normal')
    return(1-pnorm(Qclust,mean=fit$estimate[1],sd=fit$estimate[2]))
  }
}

## Test significance of a split in a hierarchical tree of pre-computed clusters (parallelized).
# data = Genes x cells UMI matrix.
# ids1 = Cell indices belonging to the left branch.
# ids2 = Cell indices belonging to the right branch.
# var.genes = Informative gene names.
# num_PCs = Number of principal components to use.
# cores = Number of cores to use.

test_split_posthoc.parallel <- function(data,ids1,ids2,var.genes,num_PCs,cores,alpha_level) { 
  # Re-order data
  cell1s <- data[,ids1]
  cell2s <- data[,ids2]
  true <- cbind(cell1s,cell2s)
  
  # Run dimension reduction on just these clusters
  pdev <- poisson.dev(true[var.genes,])
  pdev <- t(scale(t(pdev),scale=F))
  gm_sub <- RSpectra::eigs_sym(as.matrix(tcrossprod(pdev)),k=num_PCs)
  gm_sub.x <- t(crossprod(gm_sub$vectors,pdev))
  gm_sub.d <- dist(gm_sub.x[,1:2])
  
  # Compute average silhouette
  labs <- c(rep(1,length(ids1)),rep(2,length(ids2)))
  s <- silhouette(as.integer(labs),dist=gm_sub.d)
  Qclust <- mean(s[,3])
  
  # Determine the set of "on" genes
  phi_stat <- poisson_dispersion_stats(true[var.genes,])
  on_genes <- which(pnorm(phi_stat,lower.tail=F)<0.05)
  
  # Compute sample moments of the on genes
  on_counts <- t(true[var.genes[on_genes],])
  cov <- cov(on_counts)
  means <- colMeans(on_counts)
  
  # Use method of moments to estimate parameters
  sigmas <- log(((diag(cov)-means)/means^2)+1)
  mus <- log(means)-0.5*sigmas
  mus.sum <- tcrossprod(array(mus,dim=c(length(mus),1)),array(1,dim=c(length(mus),1)))+
    tcrossprod(array(1,dim=c(length(mus),1)),array(mus,dim=c(length(mus),1)))
  sigmas.sum <- tcrossprod(array(sigmas,dim=c(length(sigmas),1)),array(1,dim=c(length(sigmas),1)))+
    tcrossprod(array(1,dim=c(length(sigmas),1)),array(sigmas,dim=c(length(sigmas),1)))
  rhos <- log(cov/(exp(mus.sum+0.5*sigmas.sum))+1)
  rhos[is.na(rhos)|is.infinite(rhos)] <- (-10)
  diag(rhos) <- sigmas
  
  on_cov_eigs <- RSpectra::eigs_sym(as.matrix(rhos),k=num_PCs)
  num_pos <- sum(on_cov_eigs$values>0)
  on_cov.sub <- on_cov_eigs$vectors[,1:num_pos]%*%sqrt(diag(on_cov_eigs$values[1:num_pos]))
  on_cov <- tcrossprod(on_cov.sub)
  diag(on_cov) <- diag(rhos)
  on_cov <- sfsmisc::posdefify(on_cov)
  on_cov.sqrt <- t(chol(on_cov))
  on_means <- mus
  lambdas <- rowMeans(true[var.genes,])
  
  nm <- min(5000,ncol(true))
  
  if (nm==5000) {
    Qclusts2_1 <- mclapply(1:10,function(i) {
      # Generate a null dataset under the model
      null <- array(0,dim=c(length(var.genes),nm))
      rownames(null) <- names(lambdas) <- var.genes
      null[-on_genes,] <- array(rpois(ncol(null)*length(var.genes[-on_genes]),lambdas[-on_genes]),
                                dim=c(length(var.genes[-on_genes]),ncol(null)))
      Y <- exp(sweep(on_cov.sqrt%*%array(rnorm(nm*length(on_genes)),
                                         dim=c(length(on_genes),nm)),1,on_means,'+'))
      null[on_genes,] <- array(rpois(length(Y),Y),dim=dim(Y))
      
      # Compute null deviance residuals, project onto PCs, and cluster
      pdev2 <- poisson.dev(null)
      pdev2 <- t(scale(t(pdev2),scale=F))
      gm2 <- t(crossprod(gm_sub$vectors,pdev2))
      orig <- gm_sub.x
      knns <- queryKNN(orig, gm2,
                       k=15, BNPARAM=AnnoyParam())$index
      neighbor.labels <- apply(knns,1,function(x) labs[x])
      hc2 <- unlist(apply(neighbor.labels,2,function(x)
        sample(names(table(x))[as.vector(table(x))==max(table(x))],1)))
      
      # Compute average silhouette
      if (length(unique(hc2))==1|sum(hc2==1)==1|sum(hc2==2)==1) {
        0
      } else {
        s <- silhouette(as.integer(hc2),dist=dist(gm2[,1:2]))
        mean(s[,3])
      }
    },mc.cores=cores)
    
    # Quit early if p-value is much smaller than alpha level
    Qclusts2 <- unlist(Qclusts2_1)
    fit <- fitdistr(Qclusts2,'normal')
    pval <- 1-pnorm(Qclust,mean=fit$estimate[1],sd=fit$estimate[2])
    if (pval < 0.1*alpha_level) {
      return(pval)
    }
    
    # Otherwise, keep going
    Qclusts2_2 <- mclapply(11:50,function(i) {
      # Generate a null dataset under the model
      null <- array(0,dim=c(length(var.genes),nm))
      rownames(null) <- names(lambdas) <- var.genes
      null[-on_genes,] <- array(rpois(ncol(null)*length(var.genes[-on_genes]),lambdas[-on_genes]),
                                dim=c(length(var.genes[-on_genes]),ncol(null)))
      Y <- exp(sweep(on_cov.sqrt%*%array(rnorm(nm*length(on_genes)),
                                         dim=c(length(on_genes),nm)),1,on_means,'+'))
      null[on_genes,] <- array(rpois(length(Y),Y),dim=dim(Y))
      
      # Compute null deviance residuals, project onto PCs, and cluster
      pdev2 <- poisson.dev(null)
      pdev2 <- t(scale(t(pdev2),scale=F))
      gm2 <- t(crossprod(gm_sub$vectors,pdev2))
      orig <- gm_sub.x
      knns <- queryKNN(orig, gm2,
                       k=15, BNPARAM=AnnoyParam())$index
      neighbor.labels <- apply(knns,1,function(x) labs[x])
      hc2 <- unlist(apply(neighbor.labels,2,function(x)
        sample(names(table(x))[as.vector(table(x))==max(table(x))],1)))
      
      # Compute average silhouette
      if (length(unique(hc2))==1|sum(hc2==1)==1|sum(hc2==2)==1) {
        0
      } else {
        s <- silhouette(as.integer(hc2),dist=dist(gm2[,1:2]))
        mean(s[,3])
      }
    },mc.cores=cores)
    
    # Compute smoothed p-value
    Qclusts2 <- c(Qclusts2,unlist(Qclusts2_2))
    fit <- fitdistr(Qclusts2,'normal')
    return(1-pnorm(Qclust,mean=fit$estimate[1],sd=fit$estimate[2]))
  } else {
    Qclusts2 <- mclapply(1:50,function(i) {
      # Generate a null dataset under the model
      null <- array(0,dim=c(length(var.genes),nm))
      rownames(null) <- names(lambdas) <- var.genes
      null[-on_genes,] <- array(rpois(ncol(null)*length(var.genes[-on_genes]),lambdas[-on_genes]),
                                dim=c(length(var.genes[-on_genes]),ncol(null)))
      Y <- exp(sweep(on_cov.sqrt%*%array(rnorm(nm*length(on_genes)),
                                         dim=c(length(on_genes),nm)),1,on_means,'+'))
      null[on_genes,] <- array(rpois(length(Y),Y),dim=dim(Y))
      
      # Compute null deviance residuals, project onto PCs, and cluster
      pdev2 <- poisson.dev(null)
      pdev2 <- t(scale(t(pdev2),scale=F))
      gm2 <- t(crossprod(gm_sub$vectors,pdev2))
      orig <- gm_sub.x
      knns <- queryKNN(orig, gm2,
                       k=15, BNPARAM=AnnoyParam())$index
      neighbor.labels <- apply(knns,1,function(x) labs[x])
      hc2 <- unlist(apply(neighbor.labels,2,function(x)
        sample(names(table(x))[as.vector(table(x))==max(table(x))],1)))
      
      # Compute average silhouette
      if (length(unique(hc2))==1|sum(hc2==1)==1|sum(hc2==2)==1) {
        0
      } else {
        s <- silhouette(as.integer(hc2),dist=dist(gm2[,1:2]))
        mean(s[,3])
      }
    },mc.cores=cores)
  }
  
  # Compute smoothed p-value
  Qclusts2 <- unlist(Qclusts2)
  fit <- fitdistr(Qclusts2,'normal')
  return(1-pnorm(Qclust,mean=fit$estimate[1],sd=fit$estimate[2]))
}

### Single-cell significance of hierarchical clustering: full clustering pipeline with built-in hypothesis testing.
## data: Genes x cells UMIs matrix. 
## alpha: Desired family-wise error rate (default 0.05).
## num_features: Number of informative genes to include (default 2500).
## num_PCs: Number of principal components to use (default 30).
## cores: Number of cores to use if parallelizing (default is detected number of cores minus one).
## parallel: Whether or not to parallelize (default is true).

scSHC <- function(data,alpha=0.05,num_features=2500,num_PCs=30,parallel=T) {
  # Get variable features
  dev <- devianceFeatureSelection(data)
  var.genes <- rownames(data)[order(dev,decreasing=T)[1:num_features]] 
  
  # Dimension reduction and clustering
  pdev <- poisson.dev(data[var.genes,])
  gm <- RSpectra::eigs_sym(as.matrix(tcrossprod(pdev)),k=num_PCs)
  gm.x <- t(crossprod(gm$vectors,pdev))
  gm.d <- dist(gm.x)
  hcl <- hclust(gm.d,method='ward.D2') 
  dend <- as.dendrogram(hcl)
  
  # Traverse the tree and store clusters and q-FWERs
  dends_to_test <- list(dend)
  clusters <- list()
  q_fwers <- list()
  while (length(dends_to_test)>0) {
    
    # Identify the two-way split
    cuts <- dendextend::cutree(dends_to_test[[1]],k=2,order_clusters_as_data=F)
    leaves <- get_leaves_attr(dends_to_test[[1]],'label')
    alpha_level <- alpha*((length(leaves)-1)/(ncol(data)-1))
    
    # Get p-value of the split
    if (sum(cuts==1)==1|sum(cuts==2)==1) {
      test <- 1
    } else {
      if (parallel) {
        test <- test_split.parallel(data,leaves[cuts==1],leaves[cuts==2],var.genes,
                                    num_PCs=min(length(cuts)-1,num_PCs),cores=7,alpha_level)
      } else {
        test <- test_split.parallel(data,leaves[cuts==1],leaves[cuts==2],var.genes,
                                    num_PCs=min(length(cuts)-1,num_PCs),cores=1,alpha_level)
      }
    }
    
    # Compare to significance threshold
    if (test < alpha_level) {
      # If significant: add left and right branches to testing stack
      left <- dends_to_test[[1]][[1]]
      right <- dends_to_test[[1]][[2]]
      dends_to_test[[length(dends_to_test)+1]] <- left
      dends_to_test[[length(dends_to_test)+1]] <- right
      
      # Compute q-FWER
      q_fwers[[length(q_fwers)+1]] <- min(round(test*(ncol(data)-1)/(length(leaves)-1),2),1)
    } else {
      # If not significant: create cluster
      clusters[[length(clusters)+1]] <- leaves
      
      # Compute q-FWER
      q_fwers[[length(q_fwers)+1]] <- paste0('Cluster ',length(clusters),': ',
                                             min(round(test*(ncol(data)-1)/(length(leaves)-1),2),1))
    }
    dends_to_test[[1]] <- NULL
  }
  
  # Produce vector of cluster labels
  cluster_labels <- rep(0,ncol(data))
  names(cluster_labels) <- colnames(data)
  for (i in 1:length(clusters)) {
    cluster_labels[clusters[[i]]] <- i
  }
  
  return(cluster_labels)
}

### Significance analysis of pre-computed clusters.
## data = Genes x cells UMIs matrix.
## cluster_ids = Vector of cluster labels in the same order as data columns.
## var.genes = (optional) List of gene names used for clustering.
## alpha = Desired family-wise error rate (default 0.05).
## num_features = Number of informative genes to include, if var.genes not provided (default 2500).
## num_PCs = Number of principal components to use (default 30).
## cores: Number of cores to use if parallelizing (default is detected number of cores minus one).
## parallel: Whether or not to parallelize (default is true).

testClusters <- function(data,cluster_ids,var.genes=NULL,
                         alpha=0.05,num_features=2500,num_PCs=30,parallel=T) {
  # Get variable features
  if (is.null(var.genes)) {
    dev <- devianceFeatureSelection(data)
    var.genes <- rownames(data)[order(dev,decreasing=T)[1:num_features]] 
  } else {
    var.genes <- intersect(var.genes,rownames(data))
  }
  
  # Create hierarchy of clusters
  pseudobulk <- aggregate(t(data[var.genes,]),by=list(cluster_ids),sum)
  rownames(pseudobulk) <- pseudobulk[,1]
  pseudobulk <- pseudobulk[,2:ncol(pseudobulk)]
  pseudobulk <- sweep(pseudobulk,1,rowSums(pseudobulk),'/')
  pseudobulk.d <- dist(pseudobulk)
  pseudobulk.hc <- hclust(pseudobulk.d,method='ward.D2')
  dend <- as.dendrogram(pseudobulk.hc)
  
  # Traverse the tree and store clusters and q-FWERs
  dends_to_test <- list(dend)
  q_fwers <- list()
  clusters <- list()
  while (length(dends_to_test)>0) {
    
    # Identify the two-way split
    cuts <- dendextend::cutree(dends_to_test[[1]],k=2,order_clusters_as_data=F)
    leaves <- get_leaves_attr(dends_to_test[[1]],'label')
    ids1 <- which(cluster_ids%in%leaves[cuts==1])
    ids2 <- which(cluster_ids%in%leaves[cuts==2])
    alpha_level <- alpha*((sum(cluster_ids%in%leaves)-1)/(ncol(data)-1))
    
    # Get p-value of the split
    if (parallel) {
      test <- test_split_posthoc.parallel(data,ids1,ids2,var.genes,
                         num_PCs=min(length(c(ids1,ids2))-1,num_PCs),cores=7,alpha_level)
    } else {
      test <- test_split_posthoc.parallel(data,ids1,ids2,var.genes,
                         num_PCs=min(length(c(ids1,ids2))-1,num_PCs),cores=1,alpha_level)
    }
    
    # Compare to significance threshold
    if (test < alpha_level) {
      # If significant: get left and right branches
      left <- dends_to_test[[1]][[1]]
      right <- dends_to_test[[1]][[2]]
      
      # If there are no further splits down a branch, make it a cluster
      # Otherwise, add it to the stack
      if (max(get_branches_heights(left))!=(-Inf)) {
        dends_to_test[[length(dends_to_test)+1]] <- left
      } else {
        clusters[[length(clusters)+1]] <- get_leaves_attr(left,'label')
      }
      if (max(get_branches_heights(right))!=(-Inf)) {
        dends_to_test[[length(dends_to_test)+1]] <- right
      } else {
        clusters[[length(clusters)+1]] <- get_leaves_attr(right,'label')
      }
      
      # Compute q-FWER
      q_fwers[[length(q_fwers)+1]] <- min(round(test*(ncol(data)-1)/(sum(cluster_ids%in%leaves)-1),2),1)
    } else {
      # If not significant: create cluster
      clusters[[length(clusters)+1]] <- leaves
      
      # Compute q-FWER
      q_fwers[[length(q_fwers)+1]] <- paste0('Cluster ',length(clusters),': ',
                                             min(round(test*(ncol(data)-1)/(sum(cluster_ids%in%leaves)-1),2),1))
    }
    dends_to_test[[1]] <- NULL
  }
  
  # Return vector of new, merged cluster labels
  cluster_labels <- cluster_ids
  for (x in 1:length(clusters)) {
    cluster_labels[which(cluster_ids%in%clusters[[x]])] <- paste0('new',x)
  }
  return(cluster_labels)
}

# Helper function: compute Poisson deviance residuals
poisson.dev <- function (y) {
  n <- colSums(y)
  pis <- rowSums(y)/sum(y)
  mu <- crossprod(array(pis,dim=c(1,length(pis))),array(n,dim=c(1,length(n))))
  d <- 2 * (y * log(ifelse(y == 0, 1, y/mu)) - (y - mu))
  d[d<0] <- 0 # sometimes numerical issues create tiny negative numbers
  sqrt(d)*ifelse(y>mu,1,-1)
}

# Helper function: compute dispersion statistics
poisson_dispersion_stats <- function(y){
  n <- colSums(y)
  pis <- rowSums(y)/sum(y)
  mu <- crossprod(array(pis,dim=c(1,length(pis))),array(n,dim=c(1,length(n))))
  y2 <- (y - mu)^2 / mu
  
  disp <- rowSums(y2)/ncol(y2)
  
  if (!'matrix'%in%class(y2)) {
    y2 <- as.matrix(y2)
  }
  
  return(sqrt(ncol(y))*(disp-1)/sqrt(matrixStats::rowVars(y2)))
}



