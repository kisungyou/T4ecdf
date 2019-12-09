#' EP-means Algorithm for Clustering Empirical Distributions
#' 
#' EP-means is a variant of k-means algorithm adapted to cluster 
#' multiple empirical cumulative distribution functions under metric structure 
#' induced by Earth Mover's Distance.
#' 
#' @param elist a length \eqn{N} list of either vector or \code{ecdf} objects.
#' @param k the number of clusters.
#' 
#' @return a named list containing \describe{
#' \item{cluster}{an integer vector indicating the cluster to which each \code{ecdf} is allocated.}
#' \item{centers}{a length \eqn{k} list of centroid \code{ecdf} objects.}
#' }
#' 
#' @examples 
#' ## 3 sets of 1d samples, 10 each and add some noise
#' #    set 1 : mixture of two gaussians
#' #    set 2 : single gamma distribution
#' #    set 3 : mixture of gaussian and gamma
#' 
#' # generate data
#' myn   = 50
#' elist = list()
#' for (i in 1:10){
#'    elist[[i]] = stats::ecdf(c(rnorm(myn, mean=-2), rnorm(myn, mean=2)))
#' }
#' for (i in 11:20){
#'    elist[[i]] = stats::ecdf(rgamma(2*myn,1))
#' }
#' for (i in 21:30){
#'    elist[[i]] = stats::ecdf(rgamma(myn,1) + rnorm(myn, mean=3))
#' }
#' 
#' 
#' # run EP-means with k clusters with k=2,3,4
#' ep2 = epmeans(elist, k=2)
#' ep3 = epmeans(elist, k=3)
#' ep4 = epmeans(elist, k=4)
#' 
#' # run EP-means with k=3 clusters 
#' epout = epmeans(elist, k=3)
#' 
#' # 2d embedding using mds
#' dmat = T4ecdf::pdist(elist, type="wasserstein", as.dist=TRUE)
#' ebd2 = stats::cmdscale(dmat, 2)
#' 
#' ## visualize
#' #  (1) show ECDF for three types of data
#' opar = par(mfrow=c(3,3))
#' plot(elist[[10]], cex=0.1, main="2 Gaussians")
#' plot(elist[[20]], cex=0.1, main="Gamma")
#' plot(elist[[30]], cex=0.1, main="Gaussian+Gamma")
#' 
#' #  (2) per-class ECDFs
#' for (k in 1:myk){
#'   idk = which(epout$cluster==k)
#'   for (i in 1:length(idk)){
#'     if (i<2){
#'       pm = paste("class ",k," (size=",length(idk),")",sep="")
#'       plot(elist[[idk[i]]], verticals=TRUE, lwd=0.25, do.points=FALSE, main=pm)
#'     } else {
#'       plot(elist[[idk[i]]], add=TRUE, verticals=TRUE, lwd=0.25, do.points=FALSE)
#'     }
#'     plot(epout$centers[[k]], add=TRUE, verticals=TRUE, lwd=2, col="red", do.points=FALSE)
#'   }
#' }
#' 
#' #  (3) 2d embedding colored class labels
#' plot(ebd2, col=ep2$cluster, main="k=2 means", pch=19)
#' plot(ebd2, col=ep3$cluster, main="k=3 means", pch=19)
#' plot(ebd2, col=ep4$cluster, main="k=4 means", pch=19)
#' par(opar)
#' 
#' @references 
#' \insertRef{henderson_ep-means:_2015}{T4ecdf}
#' 
#' @export
epmeans <- function(elist, k=2){
  ###############################################
  # Preprocessing
  clist = elist_epmeans(elist)  # will use quantized ones only / flist = elist_fform(qlist)
  myk   = round(k)
  myn   = length(clist)
  
  # Quantization
  mylength = 1000
  qseq = seq(from=1e-6, to=1-(1e-6), length.out=mylength)
  qmat = array(0,c(myn,mylength))
  for (n in 1:myn){
    qmat[n,] = as.vector(stats::quantile(clist[[n]], qseq))
  }
  
  ###############################################
  # Rcpp k-means
  tmpcpp = cpp_kmeans(qmat, myk)$means
  
  ###############################################
  # Pairwise Distance Computation
  #   wrap
  mylist1 = list()
  mylist2 = list()
  for (n in 1:myn){
    mylist1[[n]] = stats::ecdf(as.vector(qmat[n,]))
  }
  for (k in 1:myk){
    mylist2[[k]] = stats::ecdf(as.vector(tmpcpp[k,]))
  }
  #   compute pairwise distance using Earth Mover's Distance
  pdistmat = dist2_wasserstein(mylist1, mylist2, 1)
  #   index
  label = base::apply(pdistmat, 1, which.min)
  
  
  ###############################################
  # Return : we want to add 'Silhouette'
  output = list()
  output$cluster = as.integer(label)
  output$centers = mylist2
  return(output)
}


# ## personal examples
# cdf0  = stats::ecdf(rnorm(100, sd=3))    # original ECDF
# qseq  = seq(from=0,to=1,length.out=1000) # quantile sequence
# quant = stats::quantile(cdf0, qseq)
# cdf1  = stats::ecdf(quant)
# 
# par(mfrow=c(1,2))
# plot(cdf0, main="Original")
# plot(cdf1, main="Recovered")