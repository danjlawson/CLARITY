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
    if(k==1){
        A=Ysvd$u[,1:(k),drop=FALSE] 
        X=diag(Ysvd$d[1:(k)],nrow=k,ncol=k)
    }else{
        A=Ysvd$u[,1:(k),drop=FALSE]
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
#' @param Ysvd (default: NULL, meaning calculate it) the SVD of Y. If you provide a different object containing locations in the list Ysvd$u, they will be used instead. IMPORTANT: Ysvd$d must be rotated to ensure that Ysvd$u = Ysvd$v
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
    class(ret)="ClarityScan"
    return(ret)
}

