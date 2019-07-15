###############################
#' @title Exact update for X in Clarity
#'
#' @description
#' This function improves the fit of
#'      Y = || A X A^T ||_2
#' by updating X using the inverse method. X and Y are assumed to be positive, and A is positive with rows summing to 1. A has dimension N by K.
#'
#' Seeking to solve 
#'             Y = A X t(A)
#'             Y = A D t(A)
#'      t(A) Y A =  t(A) A D t(A) A
#'  G t(A) Y A G = G t(A) A D t(A) A G
#'               = D
#' where G = inverse t(A) A
#' 
#' @param A A K by K non-negative matrix of latent weights
#' @param Y An N by N non-negative matrix of similarities
#' 
#' @keywords mixture
#' @return A non-negative matrix of dimension K by K
#' @seealso c_updateX
#' @export
#' 
c_updateX_solve<-function(A,Y){
    tAAinv=try(
        solve(t(A)%*% A),silent=TRUE)
    if(class(tAAinv)=="try-error") return(tAAinv)
    Xsol=tAAinv %*% t(A) %*% Y %*% A %*% tAAinv
    Xsol
}
###############################
#' @title Robust exact update for X in Clarity
#'
#' @description
#' This function improves the fit of
#'      Y = || A X A^T ||_2
#' by updating X using the inverse (or generalised inverse) method. X and Y are assumed to be positive, and A is positive with rows summing to 1. A has dimension N by K.
#'
#' Seeking to solve 
#'             Y = A X t(A)
#'             Y = A D t(A)
#'      t(A) Y A =  t(A) A D t(A) A
#'  G t(A) Y A G = G t(A) A D t(A) A G
#'               = D
#' where G = inverse t(A) A
#' 
#' @param A A K by K non-negative matrix of latent weights
#' @param Y An N by N non-negative matrix of similarities
#' 
#' @keywords mixture
#' @return A non-negative matrix of dimension K by K
#' @seealso c_updateX
#' @export
#' 
c_Robust_updateX_solve <- function(A,Y){
    ## Robustly find a solution for X in Y =A X t(A), falling back to the generalised inverse if needed
    Xsol=try(c_updateX_solve(A,Y),silent=TRUE)
    if(class(Xsol)=="try-error"){
        tAAinv = try(MASS::ginv(t(A) %*% A), silent = TRUE)
        Xsol = tAAinv %*% t(A) %*% Y %*% A %*% tAAinv
    }
    Xsol
}

###############################
#' @title Obtain normalised calX from X
#'
#' @description
#' This function moves between a fit represented as
#'      Y = || A X A^T ||_2
#' to one represented as
#'      Y = || A calX Lambda^T ||_2
#' where calX now is on the same scale as Y, and Lambda is A with columns summing to 1.  calX is the "naturally scaled" version of X
#' 
#' We have X = calX t(C)
#' So calX = X (t(C))^{-1}
#'
#' Where t(C) is defined by t(Lambda) = t(C) t(A), the matrix which makes Lambda be A with columns summing to 1
#'
#' @param A An N by K non-negative matrix of mixture proportions
#' @param X A K by K non-negative matrix of latent weights
#' 
#' @keywords mixture
#' @return A list containing the following:
#' \itemize{
#' \item C, the matrix C above
#' \item tCinv the matrix t(C)^{-1}
#' \item calX non-negative matrix of dimension K by K
#' }
#' 
#' @export
#' 
c_calXfromX<-function(A,X){
    tcs=colSums(A)
    C=diag(1/tcs,nrow=dim(X)[1],ncol=dim(X)[2])
    tCinv=diag(tcs,nrow=dim(X)[1],ncol=dim(X)[2])
    list(C=C,tCinv=tCinv, calX=X %*% tCinv)
}

###############################
#' @title Default objective function for Clarity
#'
#' @description
#' This function evaluates
#'      || Y - A X A^T ||_2
#' 
#' @param A An N by K non-negative matrix of mixture proportions
#' @param X A K by K non-negative matrix of latent weights
#' @param Y An N by N non-negative matrix of similarities
#' 
#' @keywords mixture
#' @return A numeric value of the objective function
#' 
#' @export
#' 
c_objectivefunction=function(A,X,Y){
    Yhat=A %*% X %*% t(A) 
    sum((Y - Yhat)^2)
}

    
###############################
#' @title Score a list of Clarity objects
#'
#' @description
#' Takes a list of Clarity results and evaluate their objective function at the current settings
#' 
#' @param clist A list of Clarity objects as returned by \code{\link{Clarity_Scan}}
#' 
#' @keywords mixture
#' @return A numeric vector of length equal to the length of clist, containing the objective function evaluations
#' 
#' @export
#' 
c_listscore=function(clist){
    if(class(clist)!="ClarityScan") stop("clist must be of class ClarityScan as returned by Clarity_Scan")
    sapply(clist$scan,function(x)c_objectivefunction(x$A,x$X,x$Y))
}


###############################
#' @title Run The full Clarity algorithm for a fixed number of mixtures
#'
#' @description
#' Learn a mixture model representation of a provided similarity matrix Y using Clarity. Either uses SVD ( \code{\link{c_SVD_fixedK}}) or Multiplicative method (\code{\link{c_multiplicative_fixedK}}); see \code{\link{Clarity_Scan}} for a high level overview.
#' 
#' @param Y An N by N non-negative matrix of similarities
#' @param k an integer giving the dimension of the fit to be used
#' @param method (Default: "SVDX") Method to use. Fuzzy matched to either "SVD","SVDX","SVDmix" or "Multiplicative"
#' @param verbose (Default TRUE) whether to print progress information.
#' @param ... futher parameters to be passed to \code{\link{c_multiplicative_fixedK}} or \code{\link{c_SVD_fixedK}}
#' 
#' @keywords mixture
#' @return A list of class "Clarity" containing the following:
#' \itemize{
#' \item Y, the provided non-negative N by N similarity matrix
#' \item A, the inferred non-negative N by K mixture matrix
#' \item X, the inferred non-negative K by K latent matrix
#' \item calX, the inferred non-negative K by K 'naturally scaled' latent matrix
#' \item C, the matrix C above relating X to calX
#' \item tCinv, the inverse of C^T, used in relating X to calX
#' \item objective, the objective function evaluated at the proposed solution
#' \item prediction, A X A^T
#' \item Yresid, the residuals at the proposed solution
#' \item k, the provided number of mixture components
#' \item ..., additional content returned from the method.
#' }
#' 
#' @seealso The methods are \code{\link{c_SVD_fixedK}}) for SVD or \code{\link{c_multiplicative_fixedK}} for Multiplicative. \code{\link{Clarity_Scan}} is the recomended interface and gives a high level overview.
#' @export
#' @examples
#' \donttest{
#' scansvd=Clarity_fixedK(dataraw,3) # Use the SVD method
#' scanmult=Clarity_fixedK(dataraw,3,method="M") # Use the Multiplicative method
#' }
Clarity_fixedK <- function(Y,k,method="SVDX",verbose=TRUE,...){
    eargs=as.list(match.call())
    method=c_argmatch(method,c("SVD","Multiplicative"))
    if("Ysvd"%in%names(eargs)) {
        Ysvd=eval(eargs[["Ysvd"]])
    }else {
        if(method!="Multiplicative"){
            if(verbose)print("Computing SVD of Y. You can skip this by providing it as Ysvd")
            Ysvd=svd(Y-rowMeans(Y))
        }else Ysvd=NULL
    }
    if(method=="SVDmix"){
        return(c_SVD_fixedK(Ysvd,Y,k,verbose,...))
    }else if(method=="SVD"){
        return(c_simpleSVD_fixedK(Ysvd,Y,k,verbose=verbose,Xtype="Sigma"))
    }else if(method=="SVDX"){
        return(c_simpleSVD_fixedK(Ysvd,Y,k,verbose=verbose,Xtype="X"))
    }else if(method=="Multiplicative"){
        return(c_multiplicative_fixedK(Y,k,verbose=verbose,...))
    }else{
        stop("Invalid method")
    }
}


###############################
#' @title Predict one Similarity from a Clarity or ClarityScan object
#'
#' @description
#' Takes each Clarity object and learns a new X for the new data Ynew provided.
#' 
#' @param Ynew The data to be predicted
#' @param clist The learned Clarity object, either of class "Clarity" or of class "ClarityScan"
#' @param Ysvd Default NULL. SVD of Y, which is needed for the method="SVD". You can provide it or it will be computed for you when using this method. NB: It should be the SVD of Y-rowMeans(Y) to remove the mean effect!
#' 
#' @keywords mixture
#' @return An object of the same class as clist provided, with updated Y, X and derived features.
#' 
#' @seealso \code{\link{Clarity_Scan}} to generate an estimate for clist to be used here, \code{\link{plot.Clarity}} and \code{\link{plot.ClarityScan}} for plotting.
#' @export
#' @examples
#' \donttest{
#' scanraw=Clarity_Scan(dataraw) # Generate an initial run
#' # Apply it to a new dataset with the same structure
#' scanrepfromraw=Clarity_Predict(datarep,clist=scanraw) 
#' # Apply it to a new dataset with slightly different structure
#' scanmixfromraw=Clarity_Predict(datamix,clist=scanraw) 
#' }
Clarity_Predict<-function(Ynew,clist,Ysvd=NULL) {
    if(class(clist)=="ClarityScan") {
        if(is.null(Ysvd)) Ysvd=svd(Ynew-rowMeans(Ynew))
        ret=list()
        ret$scan=lapply(clist$scan,Clarity_Predict,Ynew=Ynew,Ysvd=Ysvd)
        ret$objectives=sapply(ret$scan,function(x)x$objective)
        ret$klist=clist$klist
        ret$Y=Ynew
        ret$Ysvd=Ysvd
        ret$kmax=clist$kmax
        class(ret)="ClarityScan"
        return(ret)
    }else if(class(clist)=="Clarity"){
        ret=clist
        if(!is.null(Ysvd)) ret$Ysvd=Ysvd # Only for book keeping, we don't need it
        ret$X=c_Robust_updateX_solve(ret$A,Ynew)
        ret$calX=c_calXfromX(ret$A,ret$X)$calX
        ret$Y=Ynew
        ret$prediction=ret$A %*% ret$X %*% t(ret$A)
        ret$Yresid=ret$Y-ret$prediction 
        ret$objective=c_objectivefunction(ret$A,ret$X,Ynew)
        return(ret)
    }else{
        stop("Must provide a Clarity or ClarityScan object")
    }
}

###############################
#' @title Extend a ClarityScan to match new data
#'
#' @description
#' Extend a ClarityScan that has been fitted to new data, to additionally fit the residuals of that data.
#' First call \code{\link{Clarity_Predict}} to fit the object to the new data, then call this to fit the residuals.
#'
#' This function is not yet ready for production.
#'
#' Currently only \code{\link{c_SVD_ResidualExtend}} is available as a method.
#' 
#' @param clist A learned Clarity object, of class "Clarity".
#' @param kmax (default: 10): The number of additional components to add
#' @param extend (default: "Residual") extension method. Currently the only one available.
#' @param method (default: "SVD") representation method. Currently the only one available.
#' @param verbose (default: TRUE) verbose output?
#' 
#' @keywords mixture
#' @return An object of class "ClarityScanExtend", which is a list containing the same items as class "ClarityScan" (see \code{\link{Clarity_Scan}}) with the following modifications:
#' \itemize{
#' \item scan: A list of length kmax, each containing an object of class "ClarityExtend" containing the same items as a "Clarity" object as returned from \item{\link{Clarity_fixedK}} with the following additional items:
#' \itemize{
#' \item k0: The initial k from which the extension was run
#' \item R: The distance matrix of the residuals, which was fitted in the extension
#' \item Rsvd: The SVD of the distance matrix of the residuals.
#' }
#' \item kmax: The provided kmax
#' \item Y: The original data to fit
#' \item R: The residual structure that was fit in the extension
#' \item Rsvd: The SVD of R
#' \item k0: the originally used k
#' \item klist: the dimensions of the objects in "scan"
#' \item objectives: the objectives of the scan
#' }
#' 
#' @seealso \code{\link{Clarity_Scan}} to generate an estimate for clist to be used here, and \code{\link{Clarity_Predict}} to make an initial prediction of the new data.
Clarity_Extend <- function(clist,kmax=10,
                           extend="Residual",method="SVD",
                           verbose=TRUE){
    method=c_argmatch(method,c("SVD"))
    extend=c_argmatch(extend,c("Residual"))
    if(method=="SVD" && extend=="Residual"){
        c_SVD_ResidualExtend(clist,kmax,verbose)
    }else{
        stop("Method/Extension not implemented")
    }
}


###############################
#' @title Run the full Clarity algorithm to fit a range of mixture models to data.
#'
#' @description
#' Learn a mixture model representation of a provided similarity or distance matrix Y using Clarity.
#'
#' Several methods are available:
#' \itemize{
#' \item "SVDX": A fast method that provides high quality predictions and zero fuss. Uses \code{\link{c_simpleSVD_Scan}} and can tolerate asymmetric Y.
#' \item "SVD": A fast method that provides high quality predictions and zero fuss. Uses \code{\link{c_simpleSVD_Scan}} using the SVD singluar value matrix for X, assuming symmetry.
#' \item "SVDmix": A relatively fast method that provides high quality predictions based on forming a mixture around an SVD. Uses \code{\link{c_SVD_Scan}}.
#' \item "Multiplicative": A slower method that tries to find a minimum volume encolosing simplex for the data, making the output into a legititate mixture solution with additional assumptions. Uses \code{\link{Clarity_Multiplicative_Scan}}.
#' }
#' 
#' @param Y The similarity or distance matrix to be fitted. In principle any square matrix is allowed but algorithmically you may have problems if it is non-negative or rank-deficient.
#' @param kmax (default: 20): The maximum number of mixture components to be considered. All values less than this are explored
#' @param method (Default: "SVD") Fuzzy matched to either "SVD" or "Multiplicative"
#' @param clist (Default: NULL, meaning start from scratch) A previously learned  "ClarityScan" object. Useful for adding additional iterations to \code{\link{Clarity_Multiplicative_Scan}}.
#' @param verbose (default: TRUE) verbose output?
#' @param ... additional parameters to the methods; principally tuning parameters for  \code{\link{Clarity_Multiplicative_Scan}}.
#' 
#' @keywords mixture
#' @return An object of class "ClarityScan", which is a list containing:
#' \itemize{
#' \item scan A list of length (kmax)-1, each containing an object of class "Clarity" as described in \code{\link{Clarity_fixedK}} (for K=1..kmax)
#' \item objectives A vector of length (kmax), the objective function at each K
#' \item klist The list of K evaluated at, K=1..kmax
#' \item Y The data Y provided
#' \item kmax The maximum K requested
#' \item ... Additional method-specific details (e.g. Ysvd for method="SVD")
#' }

#' @seealso  \code{\link{Clarity_fixedK}} for the details of what a "Clarity" object is. \code{\link{Clarity_Predict}} to make an initial prediction of the new data. \code{\link{plot.ClarityScan}} for plotting, \code{\link{Clarity_Bootstrap}} for quantifying uncertainty. \code{\link{Clarity_ObjectivePlot}} for comparing quality of fit as a function of \code{k}.
#' @export
#' @examples
#' \donttest{
#' scan=Clarity_Scan(dataraw) # SVD method model fit
#' ## Onward analysis:
#' predmix=Clarity_Predict(datamix,scan) ## Core prediction
#' ## Bootstrap based on arbitrary resampling scheme
#' scanbootstrap=Clarity_Bootstrap(scan,target=datamix,D=dataraw)
#' ## Core persistence plot
#' plot(predmix,signif=scanbootstrap)
#' }
#' \dontrun{
#' scanmult=Clarity_Scan(dataraw,method="M") # Multiplicative method model fit
#' scanmult=Clarity_Scan(dataraw,method="M",clist=scanmult) # Iterate the Multiplicative method longer
#' }
Clarity_Scan <- function(Y, kmax =20,
                         method="SVDX",
                         clist=NULL, verbose=TRUE,
                         ...){
    if(!is.null(clist)) {
        method=clist$method
    }else{
        method=c_argmatch(method,c("SVD","SVDX","SVDmix","Multiplicative"))
    }
    if(method=="SVDmix"){
        clist=c_SVD_Scan(Y,kmax,clist,verbose=verbose,...)
    }else if(method=="SVD"){
        clist=c_simpleSVD_Scan(Y,kmax,clist,verbose=verbose,Xtype="Sigma",...)
    }else if(method=="SVDX"){
        clist=c_simpleSVD_Scan(Y,kmax,clist,verbose=verbose,Xtype="X",...)
    }else if(method=="Multiplicative"){
        clist=Clarity_Multiplicative_Scan(Y,kmax,clist,verbose=verbose,...)
    }else{
        stop("Invalid method")
    }
}

