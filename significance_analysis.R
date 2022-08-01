#library(glmpca)
library(scry)
#library(FastKNN)
library(AER)
#library(lqmm)
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

test_split.parallel <- function(data,ids1,ids2,var.genes,num_PCs,cores) { 
  # Re-order data
  cell1s <- data[,ids1]
  cell2s <- data[,ids2]
  true <- cbind(cell1s,cell2s)
  
  # Re-run dimension reduction on just these clusters
  vst_sub <- sctransform::vst(true[var.genes,],method='glmGamPoi',residual_type='deviance',n_genes=NULL,verbosity=0)
  if (sum(is.na(vst_sub$y))>0) {
    vst_sub$y <- vst_sub$y[rowSums(is.na(vst_sub$y))==0,]
  }
  gm_sub <- prcomp_irlba(t(vst_sub$y),num_PCs)
  gm_sub.d <- as.matrix(dist(gm_sub$x))
  rownames(gm_sub.d) <- colnames(gm_sub.d) <- colnames(true) 
  
  # Compute average silhouette
  labs <- c(rep(1,length(ids1)),rep(2,length(ids2)))
  s <- cluster:::silhouette.default.R(as.integer(labs),dmatrix=gm_sub.d)
  Qclust <- mean(s[,3])
  
  # Determine the set of "on" genes
  on_genes <- NULL
  for (gene in 1:length(var.genes)) {
    if (sum(true[var.genes[gene],])!=0) {
      pois_glm <- glm(true[var.genes[gene],]~1,family=poisson)
      od <- dispersiontest(pois_glm,trafo=1)
      if (od$p.value<0.05) {
        on_genes <- c(on_genes,gene)
      }
    }
  }
  
  # Compute sample moments of the on genes
  on_counts <- t(true[var.genes[on_genes],])
  cov <- cov(as.matrix(on_counts))
  means <- colMeans(on_counts)
  
  # Use method of moments to estimate parameters
  sigmas <- log(((diag(cov)-means)/means^2)+1)
  mus <- log(means)-0.5*sigmas
  rhos <- array(0,dim=c(length(on_genes),length(on_genes)))
  for (x in 1:(ncol(rhos)-1)) {
    for (y in (x+1):ncol(rhos)) {
      rhos[x,y] <- log((cov[x,y])/(exp(mus[x]+mus[y]+0.5*(sigmas[x]+sigmas[y])))+1)
      rhos[y,x] <- rhos[x,y]
    }
  }
  rhos[is.na(rhos)|is.infinite(rhos)] <- (-10)
  diag(rhos) <- sigmas
  
  on_cov <- as.matrix(nearPD(rhos,keepDiag=T)$mat)
  on_means <- mus
  lambdas <- rowMeans(true[var.genes,])
  
  # Generate null distribution of test statistics
  Qclusts2 <- mclapply(1:100,function(i) {
    # Generate a null dataset under the model
    null <- array(0,dim=c(length(var.genes),ncol(true)))
    colnames(null) <- colnames(true)
    rownames(null) <- names(lambdas) <- var.genes
    for (gene in var.genes[-on_genes]) {
      null[gene,] <- rpois(ncol(null),lambdas[gene])
    }
    
    Y <- t(exp(mvrnorm(n=ncol(true),on_means,on_cov)))
    null[on_genes,] <- array(rpois(length(Y),Y),dim=dim(Y))
    
    # Reduce dimensionality of null dataset and cluster
    vst2 <- sctransform::vst(null[var.genes,],method='glmGamPoi',residual_type='deviance',n_genes=NULL,verbosity=0)
    if (sum(is.na(vst2$y))>0) {
      vst2$y <- vst2$y[rowSums(is.na(vst2$y))==0,]
    }
    gm2 <- prcomp_irlba(t(vst2$y),num_PCs)
    gm2.d <- as.matrix(dist(gm2$x))
    hc2 <- cutree(hclust(as.dist(gm2.d),method='ward.D'),2)
    
    # Compute average silhouette
    s <- cluster:::silhouette.default.R(as.integer(hc2),dmatrix=gm2.d)
    mean(s[,3])
  },mc.cores=cores)
  
  # Compute smoothed p-value
  Qclusts2 <- unlist(Qclusts2)
  fit <- fitdistr(Qclusts2,'normal')
  return(1-pnorm(Qclust,mean=fit$estimate[1],sd=fit$estimate[2]))
}

## Test significance of a split in a hierarchical tree as part of sc-SHC.
# data = Genes x cells UMI matrix.
# ids1 = Cell indices belonging to the left branch.
# ids2 = Cell indices belonging to the right branch.
# var.genes = Informative gene names.
# num_PCs = Number of principal components to use.

test_split <- function(data,ids1,ids2,var.genes,num_PCs) { 
  # Re-order data
  cell1s <- data[,ids1]
  cell2s <- data[,ids2]
  true <- cbind(cell1s,cell2s)
  
  # Re-run dimension reduction on just these clusters
  vst_sub <- sctransform::vst(true[var.genes,],method='glmGamPoi',residual_type='deviance',n_genes=NULL,verbosity=0)
  if (sum(is.na(vst_sub$y))>0) {
    vst_sub$y <- vst_sub$y[rowSums(is.na(vst_sub$y))==0,]
  }
  gm_sub <- prcomp_irlba(t(vst_sub$y),num_PCs)
  gm_sub.d <- as.matrix(dist(gm_sub$x))
  rownames(gm_sub.d) <- colnames(gm_sub.d) <- colnames(true) 
  
  # Compute average silhouette
  labs <- c(rep(1,length(ids1)),rep(2,length(ids2)))
  s <- cluster:::silhouette.default.R(as.integer(labs),dmatrix=gm_sub.d)
  Qclust <- mean(s[,3])
  
  # Determine the set of "on" genes
  on_genes <- NULL
  for (gene in 1:length(var.genes)) {
    if (sum(true[var.genes[gene],])!=0) {
      pois_glm <- glm(true[var.genes[gene],]~1,family=poisson)
      od <- dispersiontest(pois_glm,trafo=1)
      if (od$p.value<0.05) {
        on_genes <- c(on_genes,gene)
      }
    }
  }
  
  # Compute sample moments of the on genes
  on_counts <- t(true[var.genes[on_genes],])
  cov <- cov(as.matrix(on_counts))
  means <- colMeans(on_counts)
  
  # Use method of moments to estimate parameters
  sigmas <- log(((diag(cov)-means)/means^2)+1)
  mus <- log(means)-0.5*sigmas
  rhos <- array(0,dim=c(length(on_genes),length(on_genes)))
  for (x in 1:(ncol(rhos)-1)) {
    for (y in (x+1):ncol(rhos)) {
      rhos[x,y] <- log((cov[x,y])/(exp(mus[x]+mus[y]+0.5*(sigmas[x]+sigmas[y])))+1)
      rhos[y,x] <- rhos[x,y]
    }
  }
  rhos[is.na(rhos)|is.infinite(rhos)] <- (-10)
  diag(rhos) <- sigmas
  
  on_cov <- as.matrix(nearPD(rhos,keepDiag=T)$mat)
  on_means <- mus
  lambdas <- rowMeans(true[var.genes,])
  
  # Generate null distribution of test statistics
  Qclusts2 <- rep(0,100)
  for (i in 1:100) {
    # Generate a null dataset under the model
    null <- array(0,dim=c(length(var.genes),ncol(true)))
    rownames(null) <- names(lambdas) <- var.genes
    for (gene in var.genes[-on_genes]) {
      null[gene,] <- rpois(ncol(null),lambdas[gene])
    }
    
    Y <- t(exp(mvrnorm(n=ncol(true),on_means,on_cov)))
    null[on_genes,] <- array(rpois(length(Y),Y),dim=dim(Y))
    colnames(null) <- colnames(true)
    
    # Reduce dimensionality of null dataset and cluster
    vst2 <- sctransform::vst(null[var.genes,],method='glmGamPoi',residual_type='deviance',n_genes=NULL,verbosity=0)
    if (sum(is.na(vst2$y))>0) {
      vst2$y <- vst2$y[rowSums(is.na(vst2$y))==0,]
    }
    gm2 <- prcomp_irlba(t(vst2$y),num_PCs)
    gm2.d <- as.matrix(dist(gm2$x))
    hc2 <- cutree(hclust(as.dist(gm2.d),method='ward.D'),2)
    
    # Compute average silhouette
    s <- cluster:::silhouette.default.R(as.integer(hc2),dmatrix=gm2.d)
    Qclust2 <- mean(s[,3])
    Qclusts2[i] <- Qclust2
  }
  
  # Compute smoothed p-value
  fit <- fitdistr(Qclusts2,'normal')
  return(1-pnorm(Qclust,mean=fit$estimate[1],sd=fit$estimate[2]))
}

## Test significance of a split in a hierarchical tree of pre-computed clusters (parallelized).
# data = Genes x cells UMI matrix.
# ids1 = Cell indices belonging to the left branch.
# ids2 = Cell indices belonging to the right branch.
# var.genes = Informative gene names.
# num_PCs = Number of principal components to use.
# cores = Number of cores to use.

test_split_posthoc.parallel <- function(data,ids1,ids2,var.genes,num_PCs,cores) { 
  # Re-order data
  cell1s <- data[,ids1]
  cell2s <- data[,ids2]
  true <- cbind(cell1s,cell2s)
  
  # Run dimension reduction on just these clusters
  vst_sub <- sctransform::vst(true[var.genes,],method='glmGamPoi',residual_type='deviance',n_genes=NULL,verbosity=0)
  if (sum(is.na(vst_sub$y))>0) {
    vst_sub$y <- vst_sub$y[rowSums(is.na(vst_sub$y))==0,]
  }
  gm_sub <- prcomp_irlba(t(vst_sub$y),num_PCs)
  gm_sub.d <- as.matrix(dist(gm_sub$x))
  rownames(gm_sub.d) <- colnames(gm_sub.d) <- colnames(true) 
  rownames(gm_sub$rotation) <- rownames(vst_sub$y)
  
  # Compute average silhouette
  labs <- c(rep(1,length(ids1)),rep(2,length(ids2)))
  s <- cluster:::silhouette.default.R(as.integer(labs),dmatrix=gm_sub.d)
  Qclust <- mean(s[,3])
  
  # Determine the set of "on" genes
  on_genes <- NULL
  for (gene in 1:length(var.genes)) {
    if (sum(true[var.genes[gene],])!=0) {
      pois_glm <- glm(true[var.genes[gene],]~1,family=poisson)
      od <- dispersiontest(pois_glm,trafo=1)
      if (od$p.value<0.05) {
        on_genes <- c(on_genes,gene)
      }
    }
  }
  
  # Compute sample moments of the on genes
  on_counts <- t(true[var.genes[on_genes],])
  cov <- cov(on_counts)
  means <- colMeans(on_counts)
  
  # Use method of moments to estimate parameters
  sigmas <- log(((diag(cov)-means)/means^2)+1)
  mus <- log(means)-0.5*sigmas
  rhos <- array(0,dim=c(length(on_genes),length(on_genes)))
  for (x in 1:(ncol(rhos)-1)) {
    for (y in (x+1):ncol(rhos)) {
      rhos[x,y] <- log((cov[x,y])/(exp(mus[x]+mus[y]+0.5*(sigmas[x]+sigmas[y])))+1)
      rhos[y,x] <- rhos[x,y]
    }
  }
  rhos[is.na(rhos)|is.infinite(rhos)] <- (-10)
  diag(rhos) <- sigmas
  
  on_cov <- as.matrix(nearPD(rhos,keepDiag=T)$mat)
  on_means <- mus
  lambdas <- rowMeans(true[var.genes,])
  
  # Generate null distribution of test statistics
  Qclusts2 <- mclapply(1:100,function(i) {
    # Generate a null dataset under the model
    null <- array(0,dim=c(length(var.genes),ncol(true)))
    rownames(null) <- names(lambdas) <- var.genes
    for (gene in var.genes[-on_genes]) {
      null[gene,] <- rpois(ncol(null),lambdas[gene])
    }
    
    Y <- t(exp(mvrnorm(n=ncol(true),on_means,on_cov)))
    null[on_genes,] <- array(rpois(length(Y),Y),dim=dim(Y))
    colnames(null) <- colnames(true)
    
    # Reduce dimensionality of null dataset, project onto PCs, and cluster with k-NN
    vst2 <- sctransform::vst(null[var.genes,],method='glmGamPoi',residual_type='deviance',n_genes=NULL,verbosity=0)
    if (sum(is.na(vst2$y))>0) {
      vst2$y <- vst2$y[rowSums(is.na(vst2$y))==0,]
    }
    gm2 <- scale(t(vst2$y[intersect(rownames(vst2$y),rownames(vst_sub$y)),]),center=T,scale=F)%*%
      gm_sub$rotation[intersect(rownames(vst2$y),rownames(vst_sub$y)),]
    orig <- scale(t(vst_sub$y[intersect(rownames(vst2$y),rownames(vst_sub$y)),]),center=T,scale=F)%*%
      gm_sub$rotation[intersect(rownames(vst2$y),rownames(vst_sub$y)),]
    knns <- queryKNN(orig, gm2,
                     k=15, BNPARAM=AnnoyParam())$index
    neighbor.labels <- apply(knns,1,function(x) labs[x])
    hc2 <- unlist(apply(neighbor.labels,2,function(x)
      sample(names(table(x))[as.vector(table(x))==max(table(x))],1)))
    
    # Compute average silhouette
    if (length(unique(hc2))==1|sum(hc2==1)==1|sum(hc2==2)==1) {
      0
    } else {
      s <- cluster:::silhouette.default.R(as.integer(hc2),dmatrix=as.matrix(dist(gm2)))
      mean(s[,3])
    }
  },mc.cores=cores)
  
  # Compute smoothed p-value
  Qclusts2 <- unlist(Qclusts2)
  fit <- fitdistr(Qclusts2,'normal')
  return(1-pnorm(Qclust,mean=fit$estimate[1],sd=fit$estimate[2]))
}

## Test significance of a split in a hierarchical tree of pre-computed clusters.
# data = Genes x cells UMI matrix.
# ids1 = Cell indices belonging to the left branch.
# ids2 = Cell indices belonging to the right branch.
# var.genes = Informative gene names.
# num_PCs = Number of principal components to use.
# cores = Number of cores to use.

test_split_posthoc <- function(data,ids1,ids2,var.genes,num_PCs) { 
  # Re-order data
  cell1s <- data[,ids1]
  cell2s <- data[,ids2]
  true <- cbind(cell1s,cell2s)
  
  # Run dimension reduction on just these clusters
  vst_sub <- sctransform::vst(true[var.genes,],method='glmGamPoi',residual_type='deviance',n_genes=NULL,verbosity=0)
  if (sum(is.na(vst_sub$y))>0) {
    vst_sub$y <- vst_sub$y[rowSums(is.na(vst_sub$y))==0,]
  }
  gm_sub <- prcomp_irlba(t(vst_sub$y),num_PCs)
  gm_sub.d <- as.matrix(dist(gm_sub$x))
  rownames(gm_sub.d) <- colnames(gm_sub.d) <- colnames(true) 
  rownames(gm_sub$rotation) <- rownames(vst_sub$y)
  
  # Compute average silhouette
  labs <- c(rep(1,length(ids1)),rep(2,length(ids2)))
  s <- cluster:::silhouette.default.R(as.integer(labs),dmatrix=gm_sub.d)
  Qclust <- mean(s[,3])
  
  # Determine the set of "on" genes
  on_genes <- NULL
  for (gene in 1:length(var.genes)) {
    if (sum(true[var.genes[gene],])!=0) {
      pois_glm <- glm(true[var.genes[gene],]~1,family=poisson)
      od <- dispersiontest(pois_glm,trafo=1)
      if (od$p.value<0.05) {
        on_genes <- c(on_genes,gene)
      }
    }
  }
  
  # Compute sample moments of the on genes
  on_counts <- t(true[var.genes[on_genes],])
  cov <- cov(on_counts)
  means <- colMeans(on_counts)
  
  # Use method of moments to estimate parameters
  sigmas <- log(((diag(cov)-means)/means^2)+1)
  mus <- log(means)-0.5*sigmas
  rhos <- array(0,dim=c(length(on_genes),length(on_genes)))
  for (x in 1:(ncol(rhos)-1)) {
    for (y in (x+1):ncol(rhos)) {
      rhos[x,y] <- log((cov[x,y])/(exp(mus[x]+mus[y]+0.5*(sigmas[x]+sigmas[y])))+1)
      rhos[y,x] <- rhos[x,y]
    }
  }
  rhos[is.na(rhos)|is.infinite(rhos)] <- (-10)
  diag(rhos) <- sigmas
  
  on_cov <- as.matrix(nearPD(rhos,keepDiag=T)$mat)
  on_means <- mus
  lambdas <- rowMeans(true[var.genes,])
  
  # Generate null distribution of test statistics
  Qclusts2 <- rep(0,100)
  for (i in 1:100) {
    # Generate a null dataset under the model
    null <- array(0,dim=c(length(var.genes),ncol(true)))
    rownames(null) <- names(lambdas) <- var.genes
    for (gene in var.genes[-on_genes]) {
      null[gene,] <- rpois(ncol(null),lambdas[gene])
    }
      
    Y <- t(exp(mvrnorm(n=ncol(true),on_means,on_cov)))
    null[on_genes,] <- array(rpois(length(Y),Y),dim=dim(Y))
    colnames(null) <- colnames(true)
    
    # Reduce dimensionality of null dataset, project onto PCs, and cluster with k-NN
    vst2 <- sctransform::vst(null[var.genes,],method='glmGamPoi',residual_type='deviance',n_genes=NULL,verbosity=0)
    if (sum(is.na(vst2$y))>0) {
      vst2$y <- vst2$y[rowSums(is.na(vst2$y))==0,]
    }
    gm2 <- scale(t(vst2$y[intersect(rownames(vst2$y),rownames(vst_sub$y)),]),center=T,scale=F)%*%
      gm_sub$rotation[intersect(rownames(vst2$y),rownames(vst_sub$y)),]
    orig <- scale(t(vst_sub$y[intersect(rownames(vst2$y),rownames(vst_sub$y)),]),center=T,scale=F)%*%
      gm_sub$rotation[intersect(rownames(vst2$y),rownames(vst_sub$y)),]
    knns <- queryKNN(orig, gm2,
                     k=15, BNPARAM=AnnoyParam())$index
    neighbor.labels <- apply(knns,1,function(x) labs[x])
    hc2 <- unlist(apply(neighbor.labels,2,function(x)
      sample(names(table(x))[as.vector(table(x))==max(table(x))],1)))
      
    # Compute average silhouette
    if (length(unique(hc2))==1|sum(hc2==1)==1|sum(hc2==2)==1) {
      Qclusts2[i] <- 0
    } else {
      s <- cluster:::silhouette.default.R(as.integer(hc2),dmatrix=as.matrix(dist(gm2)))
      Qclust2 <- mean(s[,3])
      Qclusts2[i] <- Qclust2
    }
  }
  
  # Compute smoothed p-value
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

scSHC <- function(data,alpha=0.05,num_features=2500,num_PCs=30,cores=NULL,parallel=T) {
  # Set cores if applicable
  if (parallel&is.null(cores)) {
    cores <- detectCores()-1
  }
  
  # Get variable features
  dev <- devianceFeatureSelection(data)
  var.genes <- rownames(data)[order(dev,decreasing=T)[1:num_features]] 
  
  # Dimension reduction and clustering
  vst <- sctransform::vst(data[var.genes,],method='glmGamPoi',residual_type='deviance',n_genes=NULL) 
  if (sum(is.na(vst$y))>0) {
    vst$y <- vst$y[rowSums(is.na(vst$y))==0,]
  }
  gm <- prcomp_irlba(t(vst$y),num_PCs)
  gm.d <- dist(gm$x)
  hcl <- hclust(gm.d,method='ward.D') 
  dend <- as.dendrogram(hcl)
  
  # Traverse the tree and store clusters and q-FWERs
  dends_to_test <- list(dend)
  clusters <- list()
  q_fwers <- list()
  while (length(dends_to_test)>0) {
    
    # Identify the two-way split
    cuts <- dendextend::cutree(dends_to_test[[1]],k=2,order_clusters_as_data=F)
    leaves <- get_leaves_attr(dends_to_test[[1]],'label')
    
    # Get p-value of the split
    if (parallel) {
      test <- test_split.parallel(data,leaves[cuts==1],leaves[cuts==2],var.genes,
                                  num_PCs=min(length(cuts)-1,num_PCs),cores=cores)
    } else {
      test <- test_split(data,leaves[cuts==1],leaves[cuts==2],var.genes,
                         num_PCs=min(length(cuts)-1,num_PCs))
    }
    
    # Compare to significance threshold
    alpha_level <- alpha*((length(leaves)-1)/(ncol(data)-1))
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
                         alpha=0.05,num_features=2500,num_PCs=30,cores=NULL,parallel=T) {
  # Set cores if applicable
  if (parallel&is.null(cores)) {
    cores <- detectCores()-1
  }
  
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
  pseudobulk.hc <- hclust(pseudobulk.d,method='ward.D')
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
    
    # Get p-value of the split
    if (parallel) {
      test <- test_split_posthoc.parallel(data,ids1,ids2,
                                 var.genes,num_PCs=min(length(c(ids1,ids2))-1,num_PCs),
                                 cores=cores)
    } else {
      test <- test_split_posthoc(data,ids1,ids2,
                                 var.genes,num_PCs=min(length(c(ids1,ids2))-1,num_PCs))
    }
    
    # Compare to significance threshold
    alpha_level <- alpha*((sum(cluster_ids%in%leaves)-1)/(ncol(data)-1))
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