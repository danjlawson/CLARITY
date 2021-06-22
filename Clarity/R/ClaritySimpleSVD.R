###############################
#' @title Compute the Clarity model by SVD embedding data for a fixed size of model k
#'
#' @description
#' Given data Y, and its regular svd Ysvd, and a choice of the number of SVD dimensions k to use.
#' 
#' @param Ysvd SVD of Y
#' @param Y matrix to be fit
#' @param k the number of SVD dimensions to retain
#' @param verbose (default TRUE) whether to output progress to the terminal.
#' @param Xtype (default "X") whether to construct X by optimal solution ("X") or by SVD Sigma ("Sigma"). They are equivelent when Y is symmetric.
#' 
#' @return  A Clarity object as described in \code{\link{Clarity_fixedK}}
#' @seealso \code{\link{Clarity_fixedK}} for the recommended interface.
#' @export
#' 
c_simpleSVD_fixedK <- function(Ysvd,Y,k,verbose=TRUE,Xtype="X"){
    ## Make a solution space from the SVD of Y
    ## k is here the number of dimensions in A that are wanted
    if(verbose) print(paste("Computing Clarity using simple SVD from k =",k))
    A=Ysvd$u[,1:(k),drop=FALSE] 
    if(k==1){
        X=diag(Ysvd$d[1:(k)],nrow=k,ncol=k)
    }else{
        if(Xtype=="X"){
            X=c_Robust_updateX_solve(A,Y)
        }else{
            X=diag(Ysvd$d[1:(k)],nrow=k,ncol=k)
        }
    }
    rownames(A)=rownames(Y)
    Yhat=A %*% X %*% t(A)
    Yresid=Y-Yhat
    objective=sum((Y - Yhat)^2)
    ret=list(A=A,X=X,prediction=Yhat,
             Y=Y,k=k,Yresid=Yresid,objective=objective)
    if(Xtype=="X"){
        ret$method="SVDX"
    }else{
        ret$method="SVD"
    }
    class(ret)="Clarity"
    return(ret)
}

###############################
#' @title Compute the Clarity model by SVD embedding data for a range of model sizes
#'
#' @description
#' Given data Y, embed it into a space of dimension 1..kmax and fit it as an SVD
#' 
#' @param Y matrix to be fit
#' @param kmax (default: NULL, meaning use the data dimension) the maximum number of mixture components to use
#' @param Ysvd (default: NULL, meaning calculate it) the SVD of Y. If you provide a different object containing locations in the list Ysvd$u, they will be used instead. IMPORTANT: Ysvd$d must be rotated to ensure that Ysvd$u = Ysvd$v, as done by \code{\link{c_svd0}}.
#' @param verbose (default TRUE) whether to output progress to the terminal.
#' @param Xtype (default "X") whether to construct X by optimal solution ("X") or by SVD Sigma ("Sigma"). They are equivelent when Y is symmetric.
#' 
#' @return A ClarityScan object as described in \code{\link{Clarity_Scan}}, with the additional elements:
#' \itemize{
#' \item Ysvd The SVD that was used in the embedding.
#' }
#' @seealso \code{\link{Clarity_Scan}} for the recommended interface.
#' @export
#' 
c_simpleSVD_Scan <- function(Y,kmax=NULL,Ysvd=NULL,verbose=TRUE,Xtype="X"){
    ## Performs a scan over K using the SVD method
    if(is.null(Ysvd)) {
        if(verbose) print(paste("Computing SVD"))
        Ysvd=c_svd0(Y)
    }
    if(is.null(kmax)) kmax=dim(Y)[2]
    klist=1:kmax
    ret=list()
    ret$scan=lapply(klist,function(k) c_simpleSVD_fixedK(Ysvd,Y,k,verbose,Xtype))
    ret$objectives=sapply(ret$scan,function(x)x$objective)
    ret$klist=klist
    ret$Y=Y
    ret$kmax=kmax
    ret$Ysvd=Ysvd
    if(Xtype=="X") {
        ret$method="SVDX"
    }else {
        ret$method="SVD"
    }
    class(ret)="ClarityScan"
    return(ret)
}


###############################
#' @title Iteratively apply an SVD to remove the effect of diagonals
#'
#' @description
#' The self-similarity or difference is often problematic, or follows a different model to the pairwise values. This function applies SVD and computes residuals on the diagonal, removes these residuals, and iterates. The resulting object will minimize Y - X X^T
#' 
#' @param Y The similarity data
#' @param K The SVD complexity to retain
#' @param Tmax (Default: 25) Maximum number of iterations to process. Typically converges quite quickly.
#' @param tol (Default: 1e-5) Tolerance of the difference in the sum of the residuals on the diagonal.
#' @param verbose (Default: FALSE) Whether to report progress
#' 
#' @return  A list containing:
#' \itemize{
#' \item Yp.svd: the SVD object of (Y-diag(D));
#' \item Yp: Yprime = Y-diag(D);
#' \item D: the estimated D;
#' \item K: the provided K;
#' \item Dhistory: a matrix of dimension (t times n) of the diagonals at each iteration i;
#' \item reshistory: a vector of the squared residuals at each iteration i;
#' \item Tmax: The final Tmax, after any tolerance is applied.
#' }
#' @export
#' 
c_IterativeSVD <- function(Y,K,Tmax=25,tol=1e-5,verbose=FALSE){
    residuals <- function(Y, Yhat) { sum((unlist(Y-Yhat))^2) }
    update.diag <- function(Y, D, K){
        Y.svd = svd(Y-diag(D))
        if(K==1){
            Yhat = Y.svd$u[,1:K,drop=FALSE] %*% matrix(Y.svd$d[1:K],nrow=1,ncol=1) %*% t(Y.svd$v)[1:K,,drop=FALSE]
        }else{
            Yhat = Y.svd$u[,1:K,drop=FALSE] %*% diag(Y.svd$d[1:K]) %*% t(Y.svd$v)[1:K,,drop=FALSE]
        }
        rfull=residuals(Y-diag(D), Yhat)
        rdiag=sum(((diag(Y)-D)-diag(Yhat))^2)
        new.D <- diag(Y)-diag(Yhat)
        return(list(D.old=D,D=new.D,res=rfull,resdiag=rdiag))
    }
    if(verbose) print(paste0("Computing SVD for K=",K))
    D=rep(0,dim(Y)[1])
    Dmat=matrix(NA,nrow=Tmax,ncol=length(D))
    resvec=rep(NA,Tmax)
    for(t in 1:Tmax){
        Dlist = update.diag(Y,D,K)
        D=Dlist$D
        resvec[t]=Dlist$res
        Dmat[t,]=D
        if((t>2)&&abs(abs(resvec[t]-resvec[t-1])<tol)) {
            Tmax=t
            Dmat=Dmat[1:t,]
            resvec=resvec[1:t]
            break;
        }
    }
    if(verbose) print(paste0("Completed SVD for K=",K," using T=",Tmax))
    return(list(Yp.svd=c_svd0(Y-diag(D)),
                Yp=Y-diag(D),
                D=D,
                K=K,
                Dhistory=Dmat,
                reshistory=resvec,
                Tmax=Tmax))
}



###############################
#' @title Compute the Clarity model by SVD embedding data for a range of model sizes
#'
#' @description
#' Given data Y, embed it into a space of dimension 1..kmax and fit it as an SVD
#' 
#' @param Y matrix to be fit
#' @param kmax (default: NULL, meaning use the data dimension) the maximum number of mixture components to use
#' @param Ysvd (default: NULL, meaning calculate it) the SVD of Y. If you provide a different object containing locations in the list Ysvd$u, they will be used instead. IMPORTANT: Ysvd$d must be rotated to ensure that Ysvd$u = Ysvd$v
#' @param verbose (default TRUE) whether to output progress to the terminal.
#' @param Xtype (default "X") whether to construct X by optimal solution ("X") or by SVD Sigma ("Sigma"). They are equivelent when Y is symmetric.
#' @param Tmax (default: 25) As passed to \code{\link{c_IterativeSVD}}.
#' @param tol (default: 1e-5) As passed to \code{\link{c_IterativeSVD}}.
#' 
#' @return A ClarityScan object as described in \code{\link{Clarity_Scan}}, with the additional elements:
#' \itemize{
#' \item Ysvd The list of SVDs that was used in the embedding.
#' }
#' @seealso \code{\link{Clarity_Scan}} for the recommended interface.
#' @export
#' 
c_nodiagSVD_Scan <- function(Y,kmax=NULL,Ysvd=NULL,verbose=TRUE,Xtype="X",Tmax=25,tol=1e-5){
    ## Performs a scan over K using the SVD method
    if(is.null(Ysvd)) {
        if(verbose) print(paste("Computing SVDs"))
        Ysvd=c_svdlist(Y,kmax,Tmax,tol,verbose)
    }
    if(is.null(kmax)) kmax=dim(Y)[2]
    klist=1:kmax
    ret=list()
    ret$scan=lapply(klist,function(k) {
        c_simpleSVD_fixedK(Ysvd[[k]]$Yp.svd,Ysvd[[k]]$Yp,k,verbose,Xtype)
    })
    ret$objectives=sapply(ret$scan,function(x)x$objective)
    ret$klist=klist
    ret$Y=Y
    ret$kmax=kmax
    ret$Ysvd=Ysvd
    if(Xtype=="X") {
        ret$method="SVDX"
    }else {
        ret$method="SVD"
    }
    class(ret)="ClarityScan"
    return(ret)
}

###############################
#' @title Compute the SVD of Y in convenient representation for CLARITY
#'
#' @description
#' Compute the SVD of Y, but unlike in core SVD, we use the convention that u[,i] and v[,i]
#' are positively correlated, and d can be be negative. This makes it more like the PCA represetation
#' and u==v if Y was symmetric.
#' 
#' @param Y matrix to be fit
#' 
#' @return A list containing u, d, and v as returned by svd.
#' @export
#' 

c_svd0 <-function(Y){
    Ysvd=svd(Y)
    Ys=sapply(1:dim(Ysvd$u)[1],
              function(i){
                  stats::cor(Ysvd$u[,i],Ysvd$v[,i]) >0
              })
    for(i in 1:length(Ys)){
        if(!Ys[i]) {
            Ysvd$u[,i] = -Ysvd$u[,i]
            Ysvd$d[i] = - Ysvd$d[i]
        }
    }
    Ysvd
}
###############################
#' @title Compute a set of diagonal-corrected SVDs
#'
#' @description
#' Given data Y, embed it into a space of dimension 1..kmax and fit it as an SVD
#' 
#' @param Y matrix to be fit
#' @param kmax (default: NULL, meaning use the data dimension) the maximum number of mixture components to use
#' @param Tmax (default: 25) As passed to \code{\link{c_IterativeSVD}}.
#' @param tol (default: 1e-5) As passed to \code{\link{c_IterativeSVD}}.
#' @param verbose (default TRUE) whether to output progress to the terminal.
#' 
#' @return A SVDlist object, which is a list of length kmax of SVDs as returned by \code{\link{c_IterativeSVD}}
#' @seealso \code{\link{Clarity_Scan}} for the recommended interface.
#' @export
#' 
c_svdlist=function(Y,kmax,Tmax=25,tol=1e-5,verbose=FALSE){
    Ysvd=lapply(1:kmax,c_IterativeSVD,Y=Y,Tmax=Tmax,tol=tol,verbose=verbose)
    class(Ysvd)="SVDlist"
    Ysvd
}
