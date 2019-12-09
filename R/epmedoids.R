#' Medoids
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
#' # run EP-medoids with k clusters with k=2,3,4
#' ep2 = epmedoids(elist, k=2)
#' ep3 = epmedoids(elist, k=3)
#' ep4 = epmedoids(elist, k=4)
#' 
#' # run EP-medoids with k=3 clusters 
#' epout = epmedoids(elist, k=3)
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
#' @export
epmedoids <- function(elist, k=2, type=c("KS","Lp","Wasserstein"), p=2){
  ###############################################
  # Preprocessing
  clist = elist_epmeans(elist, fname="epmedoids")
  myk   = round(k)
  
  mytype = tolower(type)
  mytype = match.arg(mytype, c("ks","lp","wasserstein"))
  myp    = as.double(p)
  
  ###############################################
  # Compute Pairwise Distance
  pmat <- T4ecdf::pdist(clist, type=mytype, p=myp, as.dist=TRUE)
  
  ###############################################
  # Run PAM Algorithm
  tmpout <- cluster::pam(pmat, k=myk)
  
  ###############################################
  # Return as of epmeans
  output = list()
  output$cluster = tmpout$clustering
  output$centers = clist[tmpout$medoids]
  return(output)
  
  return(output)
}