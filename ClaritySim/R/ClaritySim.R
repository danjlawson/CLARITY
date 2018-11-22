
###############################
#' @title simulate similarities under a coalescent model
#'
#' @description
#' Simulate a coalescent relationship between K latent objects. Then simulate N samples of those objects as dirichlet-distriubted mixtures. Finally, generate a latent structure (a K by K matrix) calX which relates the latent objects to one another using cophenetic distance on the tree, and map this into the data to create the N by N similarity matrix Y (with noise)
#' 
#' @param N Number of samples
#' @param K Dimension of the latent structure
#' @param alpha (default rep(1.5,K) ) structure of calX
#' @param sigma0 (default 0.1) Noise in the vector (Can be a vector or a scalar)
#' @param A (defaul NULL) Provide the matrix A, instead of simulating it. If NULL, it is generated from rdirichlet(N,alpha) if alpha>0 
#' 
#' @keywords mixture
#' @return A list containing:
#' \itemize{
#' \item Y, an N by N non-negative matrix of similarities
#' \item Y0, Y but without the noise added
#' \item A, an N by K non-negative matrix of mixtures
#' \item X, a K by K non-negative matrix of similarities between latent objects
#' \item calX, a K by K non-negative matrix of similarities between naturally scaled latent objects
#' \item tree, the tree that was simulated in \code{\link{rcoal}} output format
#' \item sigma0 as input
#' \item alpha as input
#' }
#' @seealso \code{\link{mixCoalescent}}, \code{\link{transformCoalescent}}
#' @export
#' @examples
#' set.seed(1)
#' mysim<-simulateCoalescent(100,10) # Simlulate 100 objects in a 10 dimensional latent space
#' # This is the data dataraw from Package Clarity
#' 
simulateCoalescent=function(N, # Number of individuals
                            K, # Dimension of the latent structure
                      alpha=rep(0.2,K), # Dirichlet hyperparameter
                      sigma0=0.01,
                      A=NULL) # Noise in the population vector (Can be a vector or a scalar)
{
    tc=ape::rcoal(K)
    td=ape::cophenetic.phylo(tc)
    if(is.null(A)){
    if(any(alpha>0)){
        A=gtools::rdirichlet(N,alpha)
        to=order(apply(A,1,function(x)which(x==max(x))))
        A=A[to,]
    }else{
        A=t(sapply(1:N,function(i){
            ret=rep(0,K)
            ret[sample(1:K,1)]=1
            ret
        }))
        to=order(apply(A,1,function(x)which(x==1)))
        A=A[to,]
    }
    
    }
    calX=td
    tcs=colSums(A)
    C=diag(1/tcs)
    X=calX %*% t(C)
    Y0=A %*% X %*% t(A)
    ##    Y = Y0 * (1 + stats::rnorm(N*N,0,sigma0))
    Y = Y0  + stats::rnorm(N*N,0,sigma0)
    Y[Y<0]=0
    colnames(Y)=rownames(Y)=rownames(A)=paste0("X",1:N)
    list(Y=Y,Y0=Y0,A=A,X=X,calX=calX,
         tree=tc,sigma0=sigma0,alpha=alpha)
}

###############################
#' @title rescale similarities of a simulated coalescent to new branch lengths
#'
#' @description
#'
#' Take a simulated tree and multiply each edge length by a random amount to generate a tree with the same topology but different structure.
#' 
#' @param sim A simulated coalescent as returned by \code{\link{simulateCoalescent}}
#' @param multmin (default 0.1) minimum multiplier for branch edge lengths
#' @param multmax (default 2) maximum multiplier for branch edge lengths
#' @param standardize (default TRUE) whether the X distances should be standardized to that of the original matrix
#' 
#' @keywords mixture
#' @return A list containing the same objects as  \code{\link{simulateCoalescent}} with updated X, Y, Y0, tree elements
#' 
#' @seealso \code{\link{simulateCoalescent}}, \code{\link{mixCoalescent}}
#' @export
#' @examples
#' set.seed(1)
#' mysim<-simulateCoalescent(100,10) # Simlulate 100 objects in a 10 dimensional latent space
#' similarsim<-transformCoalescent(mysim)
#' # similarsim$Y is the data datarep from package Clarity
#' 
transformCoalescent<-function(sim,multmin=0.1,multmax=2,standardize=TRUE){
    test2=sim
    test2$tree$edge.length=test2$tree$edge.length*stats::runif(length(test2$tree$edge.length),multmin,multmax) 
    td=ape::cophenetic.phylo(test2$tree)
    A=test2$A
    X=td
    if(standardize)X=X*mean(sim$X)/mean(X)

    Y0=A %*% X %*% t(A)
    N=dim(A)[1]
#    Y = Y0 * (1 + stats::rnorm(N*N,0,test2$sigma0))
    Y = Y0  + stats::rnorm(N*N,0,test2$sigma0)
    Y[Y<0]=0
    test2$X=X
    test2$Y=Y
    test2$Y0=Y0
    test2
}


###############################
#' @title Add a mixture edge to a coalescent for use in Clarity
#'
#' @description
#'
#' Take a simulated tree as returned by \code{\link{simulateCoalescent}}, rescale its edges with \code{\link{transformCoalescent}}. Then choose a random edge, and choose another random edge at distance at least the qmin-th quantile away. Mix these edges with amount beta (so beta=0 implies not doing anything) by directly operating on A.
#' 
#' @param sim A simulated coalescent as returned by \code{\link{simulateCoalescent}}
#' @param beta (default 0.5) amount of mixture to add
#' @param qmin (default 0.5) quantile of which edges can be linked to. qmin=0 implies any edge can be chosen.
#' @param transform (default TRUE) whether to pass sim to \code{\link{transformCoalescent}} before mixing
#' @param ... extra parameters for \code{\link{transformCoalescent}}
#' 
#' @keywords mixture
#' @return A list containing the same objects as  \code{\link{simulateCoalescent}} with updated A, X, Y, Y0, tree elements
#' 
#' @seealso \code{\link{simulateCoalescent}}, \code{\link{transformCoalescent}}
#' @export
#' @examples
#' set.seed(1)
#' mysim<-simulateCoalescent(100,10) # Simlulate 100 objects in a 10 dimensional latent space
#' similarsim<-transformCoalescent(mysim)
#' alternatesim<-mixCoalescent(mysim)
#' # alternatesim$Y is the data datamix from package Clarity
#' 

mixCoalescent<-function(sim, beta=0.5,qmin=0.5,transform=TRUE,...){
    N=dim(sim$A)[1]
    K=dim(sim$A)[2]
    if(transform)test2=transformCoalescent(sim,...)
    else test2=sim
    ##    tmix=simulateCoalescent(N,K,sim$alpha,sim$sigma0)
    X=test2$X
    A=test2$A
    testi=sample(1:K,1)
    talt=which(test2$X[testi,]>=stats::quantile(test2$X[testi,],qmin))
    testi=c(testi,sample(talt,1))
    testix=which(A[,testi[1]]>0)
    A[testix,testi[2]]=(beta) * A[testix,testi[1]]
    A[testix,testi[1]]=A[testix,testi[1]] * (1-beta)
##    X=beta * tmix$X + (1-beta)*sim$X
    Y0=A %*% X %*% t(A)
    N=dim(A)[1]
#    Y = Y0 * (1 + stats::rnorm(N*N,0,test2$sigma0))
    Y = Y0  + stats::rnorm(N*N,0,test2$sigma0)
    Y[Y<0]=0
    test2$A=A
    test2$X=X
    test2$Y=Y
    test2$Y0=Y0
#    test2$tree2=tmix$tree
    test2$beta=beta
    test2
}

