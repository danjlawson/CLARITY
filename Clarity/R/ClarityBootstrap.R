###############################
#' @title Compute the standardized divergence between a prediction and the model from which it was built
#'
#' @description
#' Compute D(Y_1,Y_2,K) - (Y_1-Y_2) in the notation of Lawson et al 2019 CLARITY paper
#' 
#' @param y2pred A ClarityScan object produced by Clarity_Predict, predicted from y1scan
#' @param y1scan A ClarityScan object produced by Clarity_Scan
#' 
#' @return A list of matrices, each containing the Standardized Divergence Matrix
#' @seealso \code{\link{c_listSymmetricDivergence}} to compute the difference betwe
#' @export
#' 
c_listDivergence=function(y2pred,y1scan){
    allD12=lapply(1:length(y1scan$scan),function(i){
        (y1scan$scan[[i]]$prediction-y2pred$scan[[i]]$prediction)-(y1scan$Y-y2pred$Y)
    })
    allD12
}
###############################
#' @title Compute the standardized symmetric divergence between a prediction and the model from which it was built
#'
#' @description
#' Compute
#' [D(Y_1,Y_2,K) - (Y_1-Y_2)]-[D(Y_2,Y_1,K) - (Y_2-Y_1)]
#' = D(Y_1,Y_2,K) - D(Y_2,Y_1,K)
#' in the notation of Lawson et al 2019 CLARITY paper
#' 
#' @param y1scan A ClarityScan object produced by Clarity_Scan
#' @param y2scan A ClarityScan object produced by Clarity_Scan
#' @param y12pred A ClarityScan object produced by Clarity_Predict, predicted from y2scan
#' @param y21pred A ClarityScan object produced by Clarity_Predict, predicted from y1scan
#' 
#' @return A list of matrices, each containing the Standardized Divergence Matrix
#' @seealso \code{\link{c_listDivergence}} which is used to compute the divergences.
#' @export
#' 
c_listSymmetricDivergence=function(y1scan,y2scan,y12pred,y21pred){
    d12=c_listDivergence(y12pred,y1scan)
    d21=c_listDivergence(y21pred,y2scan)
    ret=list()
    for(i in 1:length(d12)) ret[[i]]=d12[[i]] - d21[[i]]
    ret
}
###############################
#' @title Extract objects from a Clarity or ClarityScan object
#'
#' @description
#' Extract objects such as the Yresid residuals from a ClarityScan object. This is useful for constructing Persistences or similar. Also provides a standardized interface to extract the absolute value of residuals.
#' 
#' @param cscan A Clarity or ClarityScan object
#' @param what (default="Yresid") the object to extract
#' @param summary (default=abs) any summary to apply to the result. abs is useful to create residuals that rank naturally.
#' @param diag (default=0) What to do with the diagonal. If NULL, nothing is done; otherwise the diagonal is set to this value. The diagonal often has inflated residuals and hence should normally be set to 0 for comparison.
#' @param k (default=NULL) if specified and cscan is a ClarityScan object, only the value for the specified k is returned, allowing residuals for a single k to be extracted.
#' 
#' 
#' @return A list containing what for each k in the ClarityScan
#' @seealso \code{\link{Clarity_Persistence}} to which you can pass
#' @export
#'
Clarity_Extract=function(cscan,what="Yresid",summary=abs,diag=0,k=NULL){
    if(class(cscan)=="Clarity") {
        ret=summary(cscan[[what]])
    }else if(class(cscan)=="ClarityScan"){
        if(!is.null(k)) {
            ret=summary(cscan$scan[[k]][[what]])
        }else{
            ret = lapply(cscan$scan,function(x) summary(x[[what]]))
        }
    }else stop(paste("Invalid class",class(cscan),"for cscan in Clarity_Extract"))
    if((class(ret)=="matrix") && (!is.null(diag))) diag(ret)=diag
    if((class(ret)=="list") && (!is.null(diag))) ret=lapply(ret,function(x){diag(x)=diag;x})
    return(ret)
}
############################################
###############################
#' @title Compute the Persistence of a ClarityScan object
#' 
#' @description
#'
#' Perform an operation on a list of matrices in order to compute a Persistence Structure. By default, it computes
#' $$P_{ik} = \sum_{j=1^}N (Y_{ij} - \hat{Y}^k_{ij})^2$$
#' in the notation of Lawson et al 2019 CLARITY paper, where $Y$ is the provided data matrix and $\hat{Y}^k$ is the predicted matrix at complexity $k$.
#' 
#' @param clist A ClarityScan object, or a list of matrices
#' @param f (Default="RowSumsSquared" which generates Persistences) Either a character string matching "Frobenius" or "RowSumsSquared", or a function to compute for each matrix, which should return a vector when operating on a matrix. "Frobenius" is the Frobenius norm of the matrix (sum(x^2)), which results in a vector being returned, whereas "RowSumsSquared" evaluates as rowSums(x^2) and means a persistence matrix is returned.
#' @param what (Default="Yresid") What to extract from the ClarityScan object. "Y" and "prediction" can be extracted.
#' @param summary (Default=abs) Summary to be applied to the extracted matrices.
#' 
#' @return $f$ evaluated at each $k$, which results in an $N$ by $k$ matrix by default (where $N$ is the number of data items and $k$ is the maximum complexity)
#' @export
#' @seealso \code{\link{Clarity_Extract}}
Clarity_Persistence=function(clist,f="RowSumsSquared",what="Yresid",summary=abs){
    if(class(clist)=="ClarityScan") return(Clarity_Persistence(Clarity_Extract(clist,what,summary)))
    if(class(f)=="character") {
        if(f=="Frobenius") {
            f=function(x)sum(x^2)
        }else if(f=="RowSumsSquared"){
            f=function(x)rowSums(x^2)
        }else{
            stop(paste("Unrecognised function:",f))
        }
    }
    if(class(f)!="function") stop("Invalid function in Clarity_Persistence")
    sapply(clist,f)
}

###############################
#' @title Bootstrap resample data to obtain distribution of predictions under a specified Clarity model
#'
#' @description
#'
#' Takes a Clarity or ClarityScan object, and some feature data, and bootstrap resamples the features to result in a similarity matrix. This is predicted by the model and summarised appropriately for comparison to data. Uses \code{\link{Clarity_Persistence}} for Persistences and \code{\link{Clarity_Predict}} directly for Residuals.
#'
#' @param D an N by L matrix containing N data items and L features. The features are resampled with replacement.
#' @param clearned either a ClarityScan or Clarity object. If a Clarity object, then residuals will be bootstrapped. If a ClarityScan object then Persistence will be bootstrapped, unless a specific k is given in which case Residuals are again calculated.
#' @param nbb Number of bootstraps (default: 100)
#' @param k If specified and clearned is a ClarityScan object, residuals for this k are bootstrapped. (default: NULL, meaning use Persistences)
#' @param distfn Distance function to apply to the bootstrapped data, which should return a similarity or dissimilarity matrix. Default: c_dist, which is as,matrix(dist(x)).
#' @param summary (default=abs) Summarisation of the bootstrapped value when extracting Residuals. Default: abs, which allows residuals to be easily compared. Use summary=I to return the raw data. See \code{\link{Clarity_Extract)).
#' @param diag (default=0) diagonal value to be applied when extracting Residuals. See \code{\link{Clarity_Extract)).
#' @param seed Seed to reset the random number generator with. Default: NULL, meaning don't reset.
#' @param verbose Whether to print to screen each iteration. Default: FALSE
#' @param ... Additional parameters to \code{\link{Clarity_Persistence)) (not required for normal usage)
#' @return A list of length nbs, each containing a matrix
#' @export
Clarity_Bootstrap<-function(D,clearned,nbs=100,k=NULL,distfn=c_dist,summary=abs,diag=0,seed=NULL,verbose=FALSE,...){
    ## If provided with k, bootstraps the residuals. Otherwise bootstrap the persistence.
    if(is.null(seed))seed=1:nbs
    if(length(seed)==1) set.seed(seed)
    bspers=lapply(1:nbs,function(rep){
        if(verbose)print(paste("Iteration",rep,"of",nbs))
        if(length(seed)==nbs) set.seed(rep)
        xbs=sample(1:dim(D)[2],dim(D)[2],replace=TRUE)
        Dbs=D[,xbs,drop=FALSE]
        Ybs=distfn(Dbs)
        if(is.null(k) && class(clearned)=="ClarityScan"){
            cbs=Clarity_Predict(Ybs,clearned)
            return(Clarity_Persistence(cbs,...))
        }else{
            if(class(clearned)=="ClarityScan") cbs=Clarity_Predict(Ybs,clearned$scan[[k]])
            if(class(clearned)=="Clarity") cbs=Clarity_Predict(Ybs,clearned)
            return(Clarity_Extract(cbs,summary=summary,diag=diag))
        }            
    })
    bspers
}
###############################
#' @title Compare Observed to Bootstrap values 
#'
#' @description
#' 
#' For each element in an observed matrix, extract out the corresponding elements in a list of bootstrapped versions of that matrix. Report on the position of the observed value relative to the bootstrapped values. By default this returns the results in the form of a p-value for a test with the alternative that the observed value is larger than the bootstraps (and the null that it is not)
#'
#' Note that when population="bestrow" a conservative test is applied in which the best B bootstrap values from the entire column is used as the null distribution. This is appropriate for testing Persistences extracted via \code{\link{Clarity_Persistence}}.
#'
#' @param testbs A list of bootstrapped matrices as returned from \code{\link{Clarity_Bootsrap}}
#' @param obs A matrix of observed values
#' @param population (default="element") whether the population being compared to is the "element" X[i,j], or "bestcol": the B best values in the "column" X[,i]
#' @param comparison (default=function(x,o){sum(o>x)/(length(x)+1)}) a function to score the observed value relative to the bootstraped values.
#' @return A matrix of the same shape as obs, containing the result of comparison applied to each element
#' @export
#' @seealso \code{\link{Clarity_Bootsrap}}

Clarity_BScompare=function(testbs,obs,population="element",
                           comparison=function(x,o){sum(c(x,o)>=o)/(length(x)+1)} ) {
    if(population=="element"){
        ret=t(sapply(1:dim(obs)[1],function(i){
            sapply(1:dim(obs)[2],function(j){
                b1=sapply(testbs,function(x)x[i,j])
                comparison(b1,obs[i,j])
            })
        }))
    }else if(population=="bestcol"){
        tobs=t(obs) # Transpose before we start, since we need to iterate over columns
        ret=(sapply(1:dim(tobs)[1],function(i){
            b1=as.numeric(sapply(testbs,function(x)x[,i]))
            b1=head(sort(b1,decreasing=TRUE),length(testbs))
            sapply(1:dim(tobs)[2],function(j){
                comparison(b1,tobs[i,j])
            })
        }))
    }else(stop("Invalid population in Clarity_BScompare"))
    return(ret)
}
