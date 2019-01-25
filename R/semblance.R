#' @title Compute Semblance For a Given Input Matrix or Data Frame
#' @description Kernel methods can operate in a high-dimensional, implicit feature space with low computational cost. Here, we present
#' a rank-based Mercer kernel to compute a pair-wise similarity metric, corresponding to informative representation of data. We tailor
#' the development of a kernel to encode our prior knowledge about the data distribution over a probability space. The
#' philosophical concept behind our construction is that objects whose feature values fall on the extreme of that feature’s probability
#' mass distribution are more similar to each other, than objects whose feature values lie closer to the mean. This idea
#' represents a fundamentally novel way of assessing similarity between two observations. Our kernel (henceforth called ’Semblance’)
#' naturally lends itself to the construction of a distance metric that emphasizes features whose values lie far away
#' from the mean of their probability distribution. Semblance relies on properties empirically determined from the data and
#' does not assume an underlying distribution. The use of feature ranks on a probability space ensures that Semblance is
#' computational efficacious, robust to outliers, and statistically stable, thus making it widely applicable algorithm for
#' pattern analysis. This R package accompanies the research article "Semblance: A Data-driven Kernel Redefines the Notion of Similarity", to appear in Science Advances.
#' @param X a matrix X with n observations and m features, whose Semblance Gram Matrix is to be computed
#' @return The resultant Gram Matrix after applying Semblance kernel to the input
#' @export
#' @examples
#' # Simulation Example when the user inputs a matrix with single-cell gene expression data
#' ngenes = 10
#' ncells = 10
#' nclust = 2
#' mu=c(100, 0) #mean in cluster 1, cluster 2 for informative genes
#' sigma=c(0.01, 1) #stdev in cluster 1, cluster 2 for informative genes
#' size.rare.clust = 0.1
#' prop.info.genes = 0.2
#' n.info.genes=round(prop.info.genes*ngenes)
#' n.clust1.cells = round(ncells*size.rare.clust)
#' mu1=c(rep(mu[1]*sigma[2], n.info.genes), rep(0, ngenes-n.info.genes))
#' mu2=c(rep(mu[2]*sigma[2], n.info.genes), rep(0, ngenes-n.info.genes))
#' sig1=c(rep(sigma[1], n.info.genes), rep(1, ngenes-n.info.genes))
#' sig2=c(rep(sigma[2], n.info.genes), rep(1, ngenes-n.info.genes))
#' X=matrix(ncol=ngenes, nrow=ncells, data=0)
#' for(i in 1:n.clust1.cells){
#'   X[i,] = rnorm(ngenes, mean=mu1, sd=sig1)
#' }
#' for(i in (n.clust1.cells+1):ncells){
#'   X[i,] = rnorm(ngenes, mean=mu2, sd=sig2)
#' }
#' #Compute kernels/distances
#' rks=ranksem(X)
#' @import fields
#' @import msos
ranksem <- function(X) {
  nobs=nrow(X)
  nfeatures=ncol(X)
  K=matrix(nrow=nobs, ncol=nobs, data=0)
  Kis=array(dim=c(nobs, nobs, nfeatures))
  for(i in 1:nfeatures) {
    Kis[,,i] = computeSemblanceOneFeature(X[,i])/nobs
    K=K+Kis[,,i]
  }
  return((list(K=K/nfeatures, Kis=Kis))$K)
}

#' Compute Gini-weighted Semblance
#' @param X a matrix X with n observations and m features, whose Semblance Gram Matrix is to be computed. While computing this Gram Matrix, each feature is weighed by the Gini index for efficient feature selection.
#' @return The resultant Gini-weighted Gram Matrix after applying Semblance kernel to the input
#' @export
#' @examples
#' # Simulation Example when the user inputs a matrix with single-cell gene expression data
#' ngenes = 10
#' ncells = 10
#' nclust = 2
#' mu=c(5, 1) #mean in cluster 1, cluster 2 for informative genes
#' sigma=c(2, 1) #stdev in cluster 1, cluster 2 for informative genes
#' size.rare.clust = 0.2
#' prop.info.genes = 0.2
#' n.info.genes=round(prop.info.genes*ngenes)
#' n.clust1.cells = round(ncells*size.rare.clust)
#' mu1=c(rep(mu[1]*sigma[2], n.info.genes), rep(0, ngenes-n.info.genes))
#' mu2=c(rep(mu[2]*sigma[2], n.info.genes), rep(0, ngenes-n.info.genes))
#' sig1=c(rep(sigma[1], n.info.genes), rep(1, ngenes-n.info.genes))
#' sig2=c(rep(sigma[2], n.info.genes), rep(1, ngenes-n.info.genes))
#' X=matrix(ncol=ngenes, nrow=ncells, data=0)
#' for(i in 1:n.clust1.cells){
#'   X[i,] = rnorm(ngenes, mean=mu1, sd=sig1)
#' }
#' for(i in (n.clust1.cells+1):ncells){
#'   X[i,] = rnorm(ngenes, mean=mu2, sd=sig2)
#' }
#' Noise <- matrix(rnorm(prod(dim(X)), mean=2, sd=0.4), nrow = 10)
#' X = X + Noise
#' #Compute kernels/distances
#' rks=ranksem_Gini(X)
#' @import DescTools
#' @import PerformanceAnalytics
ranksem_Gini <- function(X){
  nobs=nrow(X)
  nfeatures=ncol(X)
  Kgi=matrix(nrow=nobs, ncol=nobs, data=0)
  Kis=array(dim=c(nobs, nobs, nfeatures))
  for(i in 1:nfeatures){
    Kis[,,i] = computeSemblanceOneFeature_Gini(X[,i])/nobs
    Kgi=Kgi+Kis[,,i]
  }
  list(Kgi=Kgi/nfeatures, Kis=Kis)
}

#' Make a matrix by repeating vector v into n columns
#' @param v a vector to be operated on
#' @param n number of columns the vector will be repeated over
#' @return a matrix with repeated columns
repCol <- function(v, n){
  return(matrix(nrow=length(v), ncol=n, data=rep(v, n), byrow=FALSE))
}

#' Make a matrix by repeating vector v into n rows
#' @param v a vector to be operated on
#' @param n number of rows the vector will be repeated over
#' @return a matrix with repeated rows
repRow <- function(v,n){
  return(matrix(nrow=n, ncol=length(v), data=rep(v, n), byrow=TRUE))
}

#' Make the upper triangular part the same as the lower triangular part.
#' @param m a matrix whose upper traingular part needs to be created using the lower traingular part
#' @return a matrix where the upper triangular part the same as the lower triangular part
makeUpperLower<-function(m){
  inds.upper= which(upper.tri(m))
  m.t = t(m)
  m.lower = m.t[which(upper.tri(m))]
  m[inds.upper]=m.lower
  return(m)
}

#' Compute semblance when there is only one feature, given as a vector x.
#' @param x a vector of observations for whom a given feature has been measured or estimated
#' @return a Semblance metric for only one feature measured for several observations
computeSemblanceOneFeature<-function(x){
  #nobs=nrow(X)
  nobs=length(x)
  r2 = rank(x, ties.method="max")
  r1 = rank(x, ties.method="min")
  d1 = repCol(r2, nobs)-repRow(r1, nobs)
  d2 = repRow(r2, nobs)-repCol(r1, nobs)
  d = pmax(d1, d2)+1
  return(nobs-d)
}

#' Compute semblance when there is only one feature, given as a vector x, but weight the feature by its Gini coefficient. Use for data with strictly positive values.
#' @param x a vector of observations for whom a given feature has been measured or estimated
#' @return a Semblance metric for only one feature measured for several observations
computeSemblanceOneFeature_Gini <-function(x){
  #nobs=nrow(X)
  nobs=length(x)
  r2 = rank(x, ties.method="max")
  r1 = rank(x, ties.method="min")
  d1 = repCol(r2, nobs)-repRow(r1, nobs)
  d2 = repRow(r2, nobs)-repCol(r1, nobs)
  d = pmax(d1, d2)+1
  sem=nobs-d
  return(sem*Gini(x))
}
