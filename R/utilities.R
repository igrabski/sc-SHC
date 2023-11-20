#' Compute Poisson deviances
#'
#' @importFrom Matrix colSums rowSums
#'
#' @noRd
#'
poisson_dev_batch <- function(y,x) {
  if (is.null(x)) {
    n <- Matrix::colSums(y)
    pis <- Matrix::rowSums(y)/sum(y)
    mu <- crossprod(array(pis,dim=c(1,length(pis))),array(n,dim=c(1,length(n))))
    d <- 2 * (y * log(ifelse(y == 0, 1, y/mu)) - (y - mu))
    d[d<0] <- 0

    return(sqrt(d)*ifelse(y>mu,1,-1))
  } else {
    y1 <- lapply(unique(x),function(i) y[,x==i,drop=F])
    n <- lapply(y1,Matrix::colSums)
    pis <- lapply(y1,function(data) Matrix::rowSums(data)/sum(data))
    mu <- lapply(1:length(y1),function(ind)
      crossprod(array(pis[[ind]],dim=c(1,length(pis[[ind]]))),
                array(n[[ind]],dim=c(1,length(n[[ind]])))))
    d <- lapply(1:length(y1),function(ind)
      2 * (y1[[ind]] * log(ifelse(y1[[ind]] == 0, 1, y1[[ind]]/mu[[ind]])) -
             (y1[[ind]] - mu[[ind]])))

    res <- array(0,dim=dim(y))
    rownames(res) <- rownames(y)
    colnames(res) <- colnames(y)
    for (ind in 1:length(y1)) {
      d[[ind]][d[[ind]]<0] <- 0
      res[,x==unique(x)[ind]] <- as.matrix(sqrt(d[[ind]])*
                                             ifelse(y1[[ind]]>mu[[ind]],1,-1))
    }

    return(res)
  }
}

#' Compute Poisson dispersion statistics
#'
#' @importFrom matrixStats rowVars
#' @importFrom Matrix colSums rowSums
#'
#' @noRd
#'
poisson_dispersion_stats <- function(y){
  n <- Matrix::colSums(y)
  pis <- Matrix::rowSums(y)/sum(y)
  mu <- crossprod(array(pis,dim=c(1,length(pis))),array(n,dim=c(1,length(n))))
  y2 <- (y - mu)^2 / mu

  disp <- Matrix::rowSums(y2)/ncol(y2)

  if (!'matrix'%in%class(y2)) {
    y2 <- as.matrix(y2)
  }

  return(sqrt(ncol(y))*(disp-1)/sqrt(matrixStats::rowVars(y2)))
}

#' Perform dimension reduction
#'
#' @importFrom Matrix t
#'
#' @noRd
#'
reduce_dimension <- function(y,x,num_PCs) {
  pdev <- poisson_dev_batch(y,x)
  pdev <- t(scale(Matrix::t(pdev),scale=F))
  PCs <- RSpectra::eigs_sym(as.matrix(tcrossprod(pdev)),k=num_PCs)
  projection <- t(crossprod(PCs$vectors,pdev))

  return(list(PCs, projection))
}

#' Compute expected sum of squares from dimension reduction scores
#'
#' @noRd
#'
compute_ess <- function(redduc) {
  sum((rowSums(sweep(redduc,2,colMeans(redduc),'-')^2)))
}

#' Compute test statistic
#'
#' @noRd
#'
ward_linkage <- function(redduc,labels) {
  ess1 <- compute_ess(redduc[labels==1,])
  ess2 <- compute_ess(redduc[labels==2,])
  ess <- compute_ess(redduc)
  return((ess-(ess1+ess2))/length(labels))
}
