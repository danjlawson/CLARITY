###############################
#' @title Standard Distance function for Clarity
#'
#' @description
#' Computes the euclidean distance between all points using as.matrix(dist(x)), enforcing the diagonal to take the minimum value of the rest of each row. 
#'
#' 
#' @param x An N by K matrix of N subjects observed at K features
#' 
#' @return An N by N matrix of the distances
#' @seealso \code{\link{resimulatedDistances}}, which computes this "correctly" in a way that is unbiased for the diagonal.
#' @export
#' 
claritysim_dist=function (x){
    ## Distance with the diagonal set to the next-lowest value
    d=as.matrix(stats::dist(x))
    for(i in 1:dim(d)[1]) d[i,i]=min(d[i,-i],d[-i,i])
    return(d)
}

###############################
#' @title Distance between a simulation and replicates of itself with the same structure
#'
#' @description
#' Computes the euclidean distance between all points by replicating the data x to x', and reporting the distance between each x and x' pair. This gives a fair representation of the diagonal.
#'
#' See \code{c_dist} in package \code{Clarity} which implements a simpler solution to the diagonal.
#' 
#' @param sim A ClaritySim object as returned by \code{simulateCoalescent}
#' @param nrep (default=100) number of replications to average over
#' @param distfn (default=as.matrix(stats::dist(x))) distance function to compute
#' 
#' @return A ClaritySim object with an updated Y
#' @seealso \code{\link{claritysim_dist}}, which does it the fast way as recommended in Clarity.
#' @examples
#' \dontrun{
#' set.seed(1)
#' mysim=simulateCoalescent(36,36,L=200,
#'                          alpha=rep(0,36),sigma0=0.1)
#' mysim2a<-resimulatedDistances(mysim,nrep=500)
#' mycols=matrix(1,nrow=dim(mysim$Y)[1],ncol=dim(mysim$Y)[2])
#' diag(mycols)=2
#' 
#' plot((mysim$Y),(mysim2a$Y),
#'      xlab="default: minimum from each row",
#'      ylab="better: simulationToDistance",col=mycols)
#' abline(a=0,b=1)
#' legend("topleft",legend=c("off-diagonal","diagonal"),
#'        text.col=1:2,bty="n")
#' }
#' @export
#' 
resimulatedDistances=function(sim,nrep=100,
                              distfn=function(x)as.matrix(stats::dist(x))){
    sim2=sim
    d=dim(sim$D)[1]
    repdata=lapply(1:nrep,function(rep){
        Drep = t(apply(sim2$A, 1, function(a) {
            stats::rnorm(sim2$L, t(a) %*% sim2$D0, sim2$sigma0)
        }))
        DT=distfn(rbind(sim$D,Drep))
        DTsubset=(DT[1:d, d+(1:d)]+
                  DT[d+(1:d),1:d] )/2
        DTsubset
    })
    Y=repdata[[1]]
    if(nrep>1) {
        for(i in 2:nrep) Y=Y+repdata[[i]]
        Y=Y/nrep
    }
    sim2$Y=Y
    sim2
}

###############################
#' @title Simulate features under a ClaritySim "mixture on a tree model"
#'
#' @description
#' Simulate traits assuming random walk using the cophenetic imposed using the tree. Specifically, data samples are generated that contain drift according to the tree, then a value created by their average mixture from that tree, plus noise
#'
#' @param sim a ClaritySim simulation containing a tree
#' @param L (default NULL meaning take from sim) number of features to simulate
#' @param sigma (default NULL meaning take from sim) the rate of drift of features per unit coalescence distance
#' @param sigma0 (default NULL meaning take from sim) the 'individual-level noise' in the trait values after accounting for individual mixtures
#'
#' @keywords mixture
#' @return A ClaritySim object, a list with updated values of
#' \itemize{
#' \item D0, the trait values associated with the tips of each branch of the tree
#' \item D, the trait values associated with each data point
#' \item Y, the matrix of distances for the samples
#' \item L, as input
#' \item sigma, as input
#' \item sigma0, as input
#' }
#' @seealso \code{\link{simulateCoalescent}}
#' @export
simData=function(sim,L=NULL,sigma=NULL,sigma0=NULL){
    if(!is.null(L))sim$L=L
    if(!is.null(sigma))sim$sigma=sigma
    if(!is.null(sigma0))sim$sigma0=sigma0
    sim$D0=sapply(1:sim$L,function(l){
        ape::rTraitCont(sim$tree,sigma=sim$sigma)
    })
    sim$D=t(apply(sim$A,1,function(a){
        stats::rnorm(sim$L,t(a) %*% sim$D0,sim$sigma0)
    }))
    sim$X=as.matrix(stats::dist(sim$D0))
    sim$Y=claritysim_dist(sim$D)
    colnames(sim$X)=rownames(sim$X)=sim$tree$tip.label
    colnames(sim$Y)=rownames(sim$Y)=rownames(sim$A)
    sim
}

###############################
#' @title simulate similarities under a coalescent model
#'
#' @description
#' Simulate a coalescent relationship between K latent objects. Then simulate N samples of those objects as dirichlet-distriubted mixtures. Finally, generate a latent structure (a K by K matrix) calX which relates the latent objects to one another using cophenetic distance on the tree, and map this into the data to create the N by N similarity matrix Y (with noise)
#' 
#' @param N Number of samples
#' @param K Dimension of the latent structure
#' @param L (default 200) Number of features to generate
#' @param alpha (default rep(1.5,K) ) structure of calX
#' @param sigma (default 0.1) the rate of drift of features per unit coalescence distance
#' @param sigma0 (default 0.01) the 'individual-level noise' in the trait values after accounting for mixture
#' @param A (default NULL) Provide the matrix A, instead of simulating it. If NULL, it is generated from rdirichlet(N,alpha) if alpha>0 
#' @param tree (default NULL) Provide the tree describing the relationship between the K latent classes
#' @param Amodel (default: "sample") either "sample" or "uniform" to determine membership of clusters when alpha=0
#' @param minedge (default: 0) Minimum branch length. Any edges shorter than minedge are set to minedge, resulting in a non-ultrametric tree.
#' 
#' @keywords mixture
#' @return A ClaritySim object, a list containing:
#' \itemize{
#' \item Y, an N by N non-negative matrix of similarities
#' \item A, an N by K non-negative matrix of mixtures
#' \item D, an N by L matrix of features
#' \item D0, a K by L matrix of features of latent objects
#' \item X, a K by K non-negative matrix of similarities between latent objects
#' \item X, a K by K non-negative matrix of tree distances between latent objects
#' \item calX, a K by K non-negative matrix of similarities between naturally scaled latent objects
#' \item tree, the tree that was simulated in \code{\link{rcoal}} output format
#' \item sigma0 as input
#' \item sigma as input
#' \item alpha as input
#' \item L as input
#' 
#' }
#' @seealso \code{\link{mixCoalescent}}, \code{\link{transformCoalescent}}
#' @export
#' @examples
#' set.seed(1)
#' # Simulate 100 objects in a 10 dimensional latent space with 200 features
#' mysim<-simulateCoalescent(100,10,200) 
#' # This is the data dataraw from Package Clarity
#' 
simulateCoalescent=function(N, # Number of individuals
                            K, # Dimension of the latent structure
                            L=200, # Number of features
                            alpha=rep(0.2,K), # Dirichlet hyperparameter
                            sigma=0.1,# Noise in the creation of features
                            sigma0=0.01,# Noise in the translation between populations and individuals
                            tree=NULL, # optional: specify a tree
                            A=NULL,# optional: specify A. must have K columns
                            Amodel="sample",
                            minedge=0)  
{
    if(is.null(tree)){
        tc=ape::rcoal(K)
        if(minedge>0)tc$edge.length=sapply(tc$edge.length,function(x)max(minedge,x))

    }else tc=tree

    td=ape::cophenetic.phylo(tc)
    if(is.null(A)) {
        if(any(alpha>0)){
            A=gtools::rdirichlet(N,alpha)
            to=order(apply(A,1,function(x)which(x==max(x))))
            A=A[to,]
        }else{
            A=t(sapply(1:N,function(i){
                ret=rep(0,K)
                if(Amodel=="sample") ret[sample(1:K,1)]=1
                else ret[1+((i-1)%%K)]=1
                ret
            }))
            to=order(apply(A,1,function(x)which(x==1)))
            A=A[to,]
        }
    }
    calX=td
    tcs=colSums(A)
    C=diag(1/tcs)
    tc$tip.label=colnames(A)=rownames(calX)=colnames(calX)=paste0("t",1:K)
    rownames(A)=paste0("X",1:N)
    sim=list(A=A,calX=calX,
             tree=tc,alpha=alpha,minedge=minedge)

    sim<-simData(sim,L,sigma,sigma0)
    class(sim)="ClaritySim"
    sim
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
#' 
#' @keywords mixture
#' @return A ClaritySim list containing the same objects as  \code{\link{simulateCoalescent}} with updated X, Y, D tree elements
#' 
#' @seealso \code{\link{simulateCoalescent}}, \code{\link{mixCoalescent}}
#' @export
#' @examples
#' set.seed(1)
#' mysim<-simulateCoalescent(100,10,200)
#' # Simlulate 100 objects in a 10 dimensional latent space with 200 features
#' similarsim<-transformCoalescent(mysim)
#' # similarsim$Y is the data datarep from package Clarity
#' 
transformCoalescent<-function(sim,multmin=0.1,multmax=2){
    sim2=sim
    sim2$tree$edge.length=sim2$tree$edge.length*stats::runif(length(sim2$tree$edge.length),multmin,multmax) 
    simData(sim2)
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
#' @param fraction (default 0.5) fraction (rounded up) of the affected cluster to affect with the mixing. When fraction=0 or 1 then technically the effect is a relationship change (although sometimes this can still only be modelled at k)
#' @param ... extra parameters for \code{\link{transformCoalescent}}
#' 
#' @keywords mixture
#' @return A list containing the same objects as  \code{\link{simulateCoalescent}} with updated A, X, Y, Y0, tree elements, and additionally:
#' \itemize{
#' dest, a list of which rows of A are `destination' nodes for the mixture
#' source, a list of which rows of A are `source' nodes for the mixture
#' edges, the named tip indexes of the (dest,source) of the mixture edge
#' beta, the provided beta
#' qmin, the provided qmin
#' }
#' 
#' @seealso \code{\link{simulateCoalescent}}, \code{\link{transformCoalescent}}
#' @export
#' @examples
#' set.seed(1)
#' mysim<-simulateCoalescent(100,10,200)
#' # Simlulate 100 objects in a 10 dimensional latent space with 200 features
#' similarsim<-transformCoalescent(mysim)
#' alternatesim<-mixCoalescent(mysim)
#' # alternatesim$Y is the data datamix from package Clarity
#' 

mixCoalescent<-function(sim, beta=0.5,qmin=0.5,transform=TRUE,fraction=0.5,...){
    if(transform)sim2=transformCoalescent(sim,...)
    else sim2=sim
    
    ## Update A
    A=sim2$A
    testi=sample(1:dim(A)[2],1)
    talt=which(sim2$X[testi,]>=stats::quantile(sim2$X[testi,],qmin))
    testi=c(testi,sample(talt,1))
    testix=which(A[,testi[1]]>0)
    if(fraction<1) testix=sort(sample(testix,size=ceiling(fraction*length(testix))))
    A[testix,testi[2]]=(beta) * A[testix,testi[1]]
    A[testix,testi[1]]=A[testix,testi[1]] * (1-beta)
    sim2$A=A

    ## Update ClaritySim object
    sim2$beta=beta
    sim2$qmin=qmin
    sim2$edges = testi
    names(sim2$edges) = sim2$tree$tip.label[sim2$edges]
    sim2$dest=which(A[,testi[1]]>0)
    sim2$source=which(A[,testi[2]]>0)
    sim2$source=sim2$source[!sim2$source%in% sim2$dest]
    simData(sim2)
}

