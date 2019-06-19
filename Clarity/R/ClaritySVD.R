##############################
#' @title Fit a simplex around a ball
#'
#' @description
#' construct regular k simplex in R^k with barycenter the origin
#' related to square of the minimum distance of barycenter to origin
#' 
#' @param k Number of dimensions of the ball
#' 
#' @return A matrix with each vertex as a column, of dimenison k+1
#' @seealso \code{\link{c_scaleSimplex}} to specify the radius
#' @export
#' 
c_regSimplex <- function(k) {
	n <- k*(k + 1) / 2
	if (k == 0) {
		B <- 0
	} else if (k == 1) {
		B <- matrix(c(1, -1), nrow = 2, byrow = TRUE)
	} else {
		A <- as.matrix(cbind(c_regSimplex(k - 1), rep(-1/sqrt(n), k)))
		B <- rbind(A, c(rep(0, k - 1), k/sqrt(n)))
	}
    return(B)
}
###############################
#' @title Fit a simplex around a ball of specified radius
#'
#' @description
#' construct regular k simplex in R^k with barycenter the origin and radius r
#' related to square of the minimum distance of barycenter to origin
#' 
#' @param k Number of dimensions of the ball
#' @param r Radius of the ball
#' 
#' @return A matrix with each vertex as a column, of dimenison k+1
#' @seealso \code{\link{c_SimplexSurrounding}} to compute the ball radius from data
#' @export
#' 
c_scaleSimplex <- function(k, r) {
	n <- k*(k + 1) / 2
	c_regSimplex(k)*r*sqrt(n)
}
###############################
#' @title Fit a simplex around data
#'
#' @description
#' construct a simplex around the data x.
#' An arbitrary orientation is chosen.
#' 
#' @param x dataset to be enclosed, with data examples in rows and features in columns
#' 
#' @return A matrix with each vertex as a column, of dimenison k+1
#' @seealso \code{\link{c_scaleSimplex}} is used to compute the simplex
#' @export
#' 
c_SimplexSurrounding <- function(x){
    xmean=colMeans(x)
    xdist=sqrt(colSums((t(x)-xmean)^2))
    t(xmean + t(c_scaleSimplex(dim(x)[2],
                             max(xdist))))
}
###############################
#' @title Compute locations of data inside a simplex
#'
#' @description
#' Given data Z and an enclosing simplex ZS, compute the locations of Z. Uses a linear model (lm) to do the computation
#' 
#' @param Z locations to be fitted, N rows with k columns
#' @param ZS enclosing simplex: k+1 rows with k columns
#' 
#' @return A mixture matrix containing N rows summing to 1, with k+1 columns
#' @export
#' 
c_GetAFromSimplex <- function(Z,ZS){
    Zmean=colMeans(Z)
    tZS=t(t(ZS)-ZS[1,])[-1,,drop=FALSE]
    tZ=t(t(Z)-ZS[1,])
    tA=t(apply(tZ,1,function(z){
        unname(stats::lm(Z~.-1,data=data.frame(Z=z,t(tZS)))$coefficients)
    }))
    if(dim(Z)[2]==1)tA=t(tA) # R coerces to the wrong thing
    Aret=cbind(1-rowSums(tA),tA)
    Aret
}
###############################
#' @title Compute the Clarity model using a mixture model surrounding an SVD embedding data for a fixed size of model k
#'
#' @description
#' Given data Y, and its mean centred svd Ysvd, and a choice of the number of mixture components k to use,.
#' 
#' @param Ysvd SVD of Y-rowMeans(Y)
#' @param Y matrix to be fit
#' @param k the number of mixture components to use
#' @param verbose (default TRUE) whether to output progress to the terminal.
#' 
#' @return  A Clarity object as described in \code{\link{Clarity_fixedK}}, with the additional elements:
#' \itemize{
#' \item Z The SVD locations of the data in k-1 dimensional space
#' \item ZS The k corners of the k-1 dimensional simplex enclosing the data
#' }
#' @seealso \code{\link{Clarity_fixedK}} for the recommended interface.
#' @export
#' 
c_SVD_fixedK <- function(Ysvd,Y,k,verbose=TRUE){
    ## Make a solution space from the SVD of Y
    ## k is here the number of dimensions in A that are wanted
    if(verbose) print(paste("Computing Clarity using SVD from k =",k))
    if(k==1) {
        Z=matrix(1,ncol=0,nrow=dim(Y)[1])
        ZS=matrix(1,ncol=1,nrow=dim(Y)[1])
        A=matrix(1,ncol=1,nrow=dim(Y)[1])
    }else{
        Z=Ysvd$u[,1:(k-1),drop=FALSE] %*% diag(sqrt(Ysvd$d[1:(k-1)]),k-1,k-1)
        ZS=c_SimplexSurrounding(Z)
        A=c_GetAFromSimplex(Z,ZS)
    }
    rownames(A)=rownames(Y)
    X=c_Robust_updateX_solve(A,Y)
    Yhat=A %*% X %*% t(A)
    Yresid=Y-Yhat
    objective=sum((Y - Yhat)^2)
    ret=list(A=A,Z=Z,ZS=ZS,X=X,prediction=Yhat,
             Y=Y,k=k,Yresid=Yresid,objective=objective,
             method="SVDmix")
    class(ret)="Clarity"
    return(ret)
}

###############################
#' @title Compute the Clarity model using a mixture model surrounding an SVD embedding data for a range of model sizes
#'
#' @description
#' Given data Y, embed it into a space of dimension 1..kmax and fit it as a mixture
#' 
#' @param Y matrix to be fit
#' @param kmax (default: NULL, meaning use the data dimension) the maximum number of mixture components to use
#' @param Ysvd (default: NULL, meaning calculate it) the SVD of Y-rowMeans(Y). If you provide a different object containing locations in the list Ysvd$u, they will be used instead.
#' @param verbose (default TRUE) whether to output progress to the terminal.
#' 
#' @return  A ClarityScan object as described in \code{\link{Clarity_Scan}}, with the additional elements:
#' \itemize{
#' \item Ysvd The SVD that was used in the embedding.
#' }
#' @seealso \code{\link{Clarity_Scan}} for the recommended interface.
#' @export
#' 
c_SVD_Scan <- function(Y,kmax=NULL,Ysvd=NULL,verbose=TRUE){
    ## Performs a scan over K using the SVD method
    if(is.null(Ysvd)) {
        if(verbose) print(paste("Computing SVD"))
        Ysvd=svd(Y-rowMeans(Y))
    }
    if(is.null(kmax)) kmax=dim(Y)[2]
    klist=1:kmax
    ret=list()
    ret$scan=lapply(klist,function(k) c_SVD_fixedK(Ysvd,Y,k,verbose))
    ret$objectives=sapply(ret$scan,function(x)x$objective)
    ret$klist=klist
    ret$Y=Y
    ret$kmax=kmax
    ret$Ysvd=Ysvd
    ret$method="SVDmix"
    class(ret)="ClarityScan"
    return(ret)
}


###############################
#' @title Extend a Clarity result by fitting to the residuals, using the SVD method and a fixed k
#'
#' @description
#' Given a Clarity object, and details of the residuals that it has failed to fit, fit an additional kextra mixture components to it to learn its unique model structure.  This uses the spectral embedding of the original fit, and the spectral embedding of the residuals, to make a feature vector that can be separated into the two components.
#'
#' @param clist the Clarity object that was learned, as return by \code{\link{Clarity_Predict}}
#' @param Rsvd The SVD of R-rowMeans(R).
#' @param R The distance matrix of the residuals, to be fit
#' @param kextra The number of additional mixture components
#' @param verbose (default TRUE) Whether to output progress information
#' 
#' @return  A ClarityExtend object as described in \code{\link{Clarity_Extend}}
#' 
#' @seealso \code{\link{Clarity_Extend}} for the recommended interface.
#' @export
#' 
c_SVD_ResidualExtendFixedK <- function(clist,Rsvd,
                                       R,kextra,verbose=TRUE){
    if(verbose) print(paste("Using Residual SVD Extension method, adding k =",kextra,"dimensions to the k =",clist$k,"provided"))
    if(kextra==0) {
        testA=clist$A
        clist$Aextra=matrix(0,nrow=dim(testA)[1],ncol=0)
        rownames(clist$Aextra)=rownames(clist$A)
        Z=clist$Z
        ZS=clist$ZS
    }else{
        Zextra=Rsvd$u[,1:(kextra),drop=FALSE] %*% diag(sqrt(Rsvd$d[1:(kextra)]),kextra,kextra)
        Z=cbind(clist$Z,Zextra)
        ZS=c_SimplexSurrounding(Z)
        A=c_GetAFromSimplex(Z,ZS)
        rownames(A)=rownames(clist$A)
        testA=A
        clist$Aextra=testA[,-(1:clist$k),drop=FALSE]
    }
    testX=c_Robust_updateX_solve(testA,clist$Y)
    clist$objective0=clist$objective
    clist$Z=Z
    clist$ZS=ZS
    clist$A=testA
    clist$X=testX
    clist$Rsvd=Rsvd
    clist$R=R
    clist$prediction=testA %*% testX %*% t(testA)
    clist$Yresid=clist$Y-clist$prediction
    clist$k=clist$k+kextra
    clist$objective=c_objectivefunction(testA,
                                        testX,clist$Y)
    class(clist)="ClarityExtend"
    clist
}

###############################
#' @title Extend a Clarity result by fitting to the residuals, using the SVD method and a range of k
#'
#' @description
#' Given a Clarity object, and details of the residuals that it has failed to fit, fit an additional 0..kmax mixture components to it to learn its unique model structure.  This uses the spectral embedding of the original fit, and the spectral embedding of the residuals, to make a feature vector that can be separated into the two components.
#'
#' @param clist the Clarity object that was learned, as return by \code{\link{Clarity_Predict}}
#' @param kmax (Default: 10) The maximum number of additional mixture components to be considered
#' @param verbose (default TRUE) Whether to output progress information
#' 
#' @return  A ClarityScanExtend object as described in \code{\link{Clarity_Extend}}
#' 
#' @seealso \code{\link{Clarity_Extend}} for the recommended interface.
#' @export
#' 
c_SVD_ResidualExtend <- function(clist,kmax=10,verbose=TRUE) {
    kmax=min(kmax,dim(clist$Y)[1] - clist$k)
    ret=list()
    ret$prediction0=clist$prediction
    ret$Y=clist$Y
    ret$residual0=ret$Y - clist$prediction
    if(verbose) print("Using Residual SVD Extension method, computing SVD of residuals")
    ret$R=as.matrix(stats::dist(ret$residual0))
    ret$Ysvd=clist$Ysvd
    ret$Rsvd=svd(ret$R-rowMeans(ret$R))
    ret$k0=clist$k
    ret$kmax=kmax
    ret$klist=clist$k + 0:kmax
    ret$scan=lapply(0:kmax,function(kextra){
        c_SVD_ResidualExtendFixedK(clist,ret$Rsvd,ret$R,kextra,verbose)
    })
    ret$objectives=sapply(ret$scan,function(x)x$objective)
    ret$methpd=clist$method
    class(ret)="ClarityScanExtend"
    ret
}
