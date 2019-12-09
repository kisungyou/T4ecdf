#' Hierarchical Agglomerative Clustering for Empirical Distributions
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
#' # run EP-hclust with 3 different methods and 2 distance metrics
#' eh2ks = ephclust(elist, method="single")
#' eh3ks = ephclust(elist, method="complete")
#' eh4ks = ephclust(elist, method="median")
#' eh2w2 = ephclust(elist, type="wasserstein", method="single")
#' eh3w2 = ephclust(elist, type="wasserstein", method="complete")
#' eh4w2 = ephclust(elist, type="wasserstein", method="median")
#' 
#' # visualize
#' opar = par(mfrow=c(2,3))
#' plot(eh2ks, main="KS+single")
#' plot(eh3ks, main="KS+complete")
#' plot(eh4ks, main="KS+median")
#' plot(eh2w2, main="2 Wasserstein+single")
#' plot(eh3w2, main="2 Wasserstein+complete")
#' plot(eh4w2, main="2 Wasserstein+median")
#' par(opar)
#' 
#' @export
ephclust <- function(elist, method=c("single","complete","average","mcquitty",
                                     "ward.D","ward.D2","centroid","median"),
                     type=c("KS","Lp","Wasserstein"), p=2){
  ###############################################
  # preprocessing
  clist    = elist_epmeans(elist, fname="ephclust")
  mymethod = match.arg(method)
  mytype = tolower(type)
  mytype = match.arg(mytype, c("ks","lp","wasserstein"))
  myp    = as.double(p)
  
  ###############################################
  # compute pairwise distance
  dmat = T4ecdf::pdist(clist, type=mytype, p=myp, as.dist=TRUE)
  
  ###############################################
  # apply 'fastcluster::hclust'
  output = fastcluster::hclust(dmat, method=mymethod)
  return(output) 
}