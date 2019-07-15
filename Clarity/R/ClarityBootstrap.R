###############################
#' @title Standard Distance function for Clarity
#'
#' @description
#' Computes the euclidean distance between all points using as.matrix(dist(x)). Provided to provide an example for user-defined dissimilarity functions.
#' 
#' @param x An N by K matrix of N subjects observed at K features
#' 
#' @return An N by N matrix of the distances
#' @seealso \code{\link{Clarity_Bootstrap}}, which takes an argument distfn which can be given this or another function that generates a similarity from feature data.
#' 
#' @export
#' 
c_dist=function(x){
    as.matrix(stats::dist(x))
}

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
#' Extract objects such as Clarity objects or residuals from a ClarityScan object. This is useful for constructing Persistences or similar. Also provides a standardized interface to extract the absolute value of residuals.
#' 
#' @param cscan A Clarity or ClarityScan object
#' @param what (default="Yresid") the object to extract. If a number and cscan is a ClarityScan object, returns the associated Clarity object. If it is the string "Rsq" return the square residuals.
#' @param summary (default=abs) any summary to apply to the result. abs is useful to create residuals that rank naturally.
#' @param diag (default=0) What to do with the diagonal. If NULL, nothing is done; otherwise the diagonal is set to this value. The diagonal often has inflated residuals and hence should normally be set to 0 for comparison.
#' @param k (default=NULL) if specified and cscan is a ClarityScan object, only the value for the specified k is returned, allowing residuals for a single k to be extracted.
#' 
#' 
#' @return A list containing what for each k in the ClarityScan
#' @seealso \code{\link{Clarity_Persistence}} to which you can pass
#' @examples
#' \donttest{
#' scan=Clarity_Scan(dataraw)
#' e1=Clarity_Extract(scan) # Get a list containing the absolute residual matrices for each k
#' e2=Clarity_Extract(scan,what=5) # get the 5th Clarity object
#' e3=Clarity_Extract(scan,summary=I,k=5) # get the 5th Residual matrix
#' e4=Clarity_Extract(scan,summary=I,k=5) # get the 5th Residual matrix
#' e5=Clarity_Extract(scan,"Rsq",k=5) # get the 5th square residual matrix
#' e6=Clarity_Extract(scan,"Rsq") # Get a list containing the R squared matrix for each k
#' }
#' @export
#' 
Clarity_Extract=function(cscan,what="Yresid",summary=abs,diag=0,k=NULL){
    if(class(cscan)=="Clarity") {
        if(what=="Rsq") {
            ret=cscan[["Yresid"]]^2
        }else ret=summary(cscan[[what]])
    }else if(class(cscan)=="ClarityScan"){
        if(class(what)=="numeric") return(cscan$scan[[what]])
        if(!is.null(k)) {
            if(what=="Rsq") {
                ret=cscan$scan[[k]][["Yresid"]]^2
            }else{
                ret=summary(cscan$scan[[k]][[what]])
            }
        }else{
            ret=lapply(cscan$scan,function(x) Clarity_Extract(x,what,summary,diag,k))
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
#' \eqn{P_{ik} = \sum_{j=1^}N (Y_{ij} - \hat{Y}^k_{ij})^2}
#' in the notation of Lawson et al 2019 CLARITY paper, where Y is the provided data matrix and \eqn{\hat{Y}^k} is the predicted matrix at complexity k.
#' 
#' @param clist A ClarityScan object, or a list of matrices
#' @param f (Default="RowSumsSquared" which generates Persistences) Either a character string matching "Frobenius" or "RowSumsSquared", or a function to compute for each matrix, which should return a vector when operating on a matrix. "Frobenius" is the Frobenius norm of the matrix (sum(x^2)), which results in a vector being returned, whereas "RowSumsSquared" evaluates as rowSums(x^2) and means a persistence matrix is returned.
#' @param what (Default="Yresid") What to extract from the ClarityScan object. "Y" and "prediction" can be extracted.
#' @param diag (Default=0) Special values to give the diagonal. NULL leaves them unchanged, 0 removes from the persisetence measure. See \code{\link{Clarity_Extract}}.
#' @param summary (Default=abs) Summary to be applied to the extracted matrices.
#' 
#' @return f evaluated at each k, which results in an N by k matrix by default (where N is the number of data items and k is the maximum complexity)
#' @examples
#' \donttest{
#' scan=Clarity_Scan(dataraw) ## Core Clarity
#' predmix=Clarity_Predict(datamix,scan)
#' p=Clarity_Persistence(scan) # extract persistence
#' p2=Clarity_Persistence(predmix) # extract persistence from ClarityScan
#' par(mfrow=c(1,2))
#' Clarity_Chart(p,main="Learned persistence") # p is a matrix 
#' Clarity_Chart(p2,main="Different topology persistence") # p2 is a matrix
#' }
#' @seealso \code{\link{Clarity_Extract}}
#' @export
Clarity_Persistence=function(clist,f="RowSumsSquared",what="Yresid",summary=abs,diag=0){
    if(class(clist)=="ClarityScan") return(Clarity_Persistence(Clarity_Extract(clist,what,summary,diag=diag)))
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
#' Takes a Clarity or ClarityScan object and create an ensemble of results for how it predicts a given target matrix. It will use a rank k representation from the Clarity object to generate the best fitting version of the target that has the learned rank k structure.
#' 
#' You can either provide a list of replicated similarity matrices, or some feature data from which bootstrap resamples can be taken from the features to result in a similarity matrix.
#'
#' The assumption is that the training data (used to created clearned) is being bootstrapped. This might have a different scale and rotation to the test data, (defined in target). Therefore, if H0="structure" then the bootstrapped data is first Procrustes transformed to best match the target using \code{\link{Clarity_Procrust}}.
#'
#' Each replicated similarity matrix is predicted by the model and summarised appropriately for comparison to data. Uses \code{\link{Clarity_Persistence}} for Persistences and \code{\link{Clarity_Predict}} directly for Residuals.
#'
#' @param clearned either a ClarityScan or Clarity object. If a Clarity object, then residuals will be bootstrapped. If a ClarityScan object then Persistence will be bootstrapped, unless a specific k is given in which case Residuals are again calculated.
#' @param target (default=NULL) The target matrix Y that will be modelled by this bootstrapping result. Must be specified unless H0="value".
#' @param D (default=NULL) An N by L matrix containing N data items and L features. If provided, bootstrapped similarity matrices Ylist are generated by resampling features with replacement.
#' @param Ylist (default=NULL) A list of N by N similarity matrices generated exactly as the original data were generated.
#' @param H0 (default="structure") The null hypothesis to be tested. The default test is that the "structure" has changed. Can also be set to "value" in which case no rotation to the target is done. (NB H1 is set in )
#' @param nbs (default=100) Number of bootstraps to use when D is provided.
#' @param k If specified and clearned is a ClarityScan object, residuals for this k are bootstrapped. (default: NULL, meaning use Persistences)
#' @param distfn Distance function to apply to the bootstrapped data, which should return a similarity or dissimilarity matrix. Default: c_dist, which is as,matrix(dist(x)).
#' @param summary (default=abs) Summarisation of the bootstrapped value when extracting Residuals. Default: abs, which allows residuals to be easily compared. Use summary=I to return the raw data. See \code{\link{Clarity_Extract}}.
#' @param diag (default=0) diagonal value to be applied when extracting Residuals. See \code{\link{Clarity_Extract}}.
#' @param seed Seed to reset the random number generator with. Default: NULL, meaning don't reset.
#' @param verbose Whether to print to screen each iteration. Default: FALSE
#' @param ... Additional parameters to \code{\link{Clarity_Persistence}} (not required for normal usage)
#' @return A list of length nbs, each containing a matrix
#' @seealso \code{\link{Clarity_Scan}} to run Clarity, \code{\link{Clarity_Extract}} to extract residuals and persistences, \code{\link{plot.ClarityScan}} and \code{\link{plot.Clarity}} for plotting.
#' @examples
#' \donttest{
#' scan=Clarity_Scan(dataraw) ## Core Clarity
#' predmix=Clarity_Predict(datamix,scan)
#' ## Simple feature-based bootstrap for persistences
#' scanbootstrap=Clarity_Bootstrap(scan,target=datamix,D=datarawD)
#'
#' ## Construct a list of externally produced resamples
#' Dreplist=lapply(1:100,function(rep){
#'     datarawD[,sample(1:dim(datarawD)[2],
#'     dim(datarawD)[2], replace = TRUE)]
#' })
#' Yreplist=lapply(Dreplist,c_dist)
#' ## NB Yreplist can have come from a complex resampling scheme
#'
#' ## Bootstrap based on arbitrary resampling scheme
#' scanbootstrap=Clarity_Bootstrap(scan,target=datamix,Ylist=Yreplist)
#'
#' ## Plot bootstrapped persistence chart
#' predplot=plot(predmix,signif=scanbootstrap)
#'
#' ## Similarly for Residuals
#' k10bootstrap=Clarity_Bootstrap(Clarity_Extract(scan,10),
#'                  target=datamix,
#'                  D=datarawD)
#' predplot=plot(Clarity_Extract(predmix,10),
#'     signif=k10bootstrap)
#' }
#' @export
Clarity_Bootstrap<-function(clearned,target=NULL,
                            D=NULL,Ylist=NULL,
                            H0="structure",
                            nbs=100,k=NULL,distfn=c_dist,
                            summary=abs,diag=0,
                            seed=NULL,verbose=FALSE,...){
    ## If provided with k, bootstraps the residuals. Otherwise bootstrap the persistence.
    if(all(is.null(D)) && all(is.null(Ylist))) stop("Must provide exactly one of D or Ylist")
    if(!all(is.null(D))){
        ## We want to use D to bootstrap the data
        if(!all(is.null(Ylist))) stop("Must provide exactly one of D or Ylist")
    }else{
        ## We want to use Ylist as provided bootstraps of the data
        nbs=length(Ylist)
    }
    if(is.null(seed))seed=1:nbs
    if(length(seed)==1) set.seed(seed)
    if(!is.null(k) && class(clearned)=="ClarityScan"){
        clearned=clearned$scan[[k]]
    }
    if(is.null(target)&&H0!="value"){
        stop("Clarity_Bootstrap requires a target matrix when H0!=\"value\". The recommended usage is to set target to be the test matrix, i.e. observed \"obs\" matrix that you will provide to Clarity_BSCompare.")
    }
    bspers=lapply(1:nbs,function(rep){
        if(verbose)print(paste("Iteration",rep,"of",nbs))
        if(all(is.null(Ylist))){
            ## Use the provided D
            if(length(seed)==nbs) set.seed(seed[rep])
            xbs=sample(1:dim(D)[2],dim(D)[2],replace=TRUE)
            Dbs=D[,xbs,drop=FALSE]
            Ybs=distfn(Dbs)
        }else Ybs=Ylist[[rep]]
        if(H0=="structure") {
            return(c_BootstrapStructure(clearned=clearned,Ybs=Ybs,target=target,
                                        summary=summary,diag=diag,transform=Clarity_Procrust,
                                        ...))
        }else if(H0=="value"){
            return(c_BootstrapValue(clearned=clearned,Ybs=Ybs,...))
        }else stop("Invalid H0! must be either \"structure\" or \"value\"")
    })
    if(class(clearned)=="Clarity") class(bspers)="ClarityBootstrap"
    if(class(clearned)=="ClarityScan") class(bspers)="ClarityScanBootstrap"
    return(bspers)
}

###############################
#' @title Compute Residuals or Persistences on bootstrapped data by scaling them to a target matrix
#'
#' @description
#' Given a resampled matrix Ybs which was originally used to make a Clarity or ClarityScan object clearned, we:
#' 1. Construct a Clarity model for Ybs 
#' 2. Use it to predict the target
#' 3. Compute the residual structure of this prediction from the original learned data
#'
#' The fit is done with a k that matches that in clearned. If clearned is of class "ClarityScan" then we fit with \code{\link{Clarity_Scan}}, otherwise we fit with \code{\link{Clarity_fixedK}}. We predict the target with \code{\link{Clarity_Predict}}.
#'
#' 
#' @param clearned either a ClarityScan or Clarity object. If a Clarity object, then residuals will be bootstrapped. If a ClarityScan object then Persistence will be bootstrapped, unless a specific k is given in which case Residuals are again calculated.
#' @param Ybs The bootstrapped matrix for which the new model will be learned
#' @param target The target matrix Y that will be modelled by this bootstrapping result.
#' @param summary (default=abs) Summarisation of the bootstrapped value when extracting Residuals. Default: abs, which allows residuals to be easily compared. Use summary=I to return the raw data. See \code{\link{Clarity_Extract}}.
#' @param diag (default=0) diagonal value to be applied when extracting Residuals. See \code{\link{Clarity_Extract}}.
#' @param transform (default=Clarity_Procrust) A function to perform a transform of the data to match the target, taking two arguments, (Ybs,target). Default: Clarity_Procrust. 
#' @param ... Additional parameters to \code{\link{Clarity_Persistence}} (not required for normal usage)
#' @return A list of residuals if clearned is of class Clarity, or a list of persistences if clearned is of class ClarityScan.
#' @seealso \code{\link{Clarity_Bootstrap}} which is the recommended way to access this code.
#' @export
c_BootstrapStructure=function(clearned,Ybs,target,summary=abs,diag=0,
                              transform=Clarity_Procrust,...){
    if(class(clearned)=="ClarityScan") {
        if(class(transform)=="function") Ybs=transform(Ybs,target)
        matlist=lapply(clearned$scan,function(x){
            c_BootstrapStructure(x,Ybs,target,summary=summary,
                                 diag=diag,transform=FALSE,...)
        })
        return(Clarity_Persistence(matlist))
    }
    ## Rank k fit to the data
    if(class(transform)=="function") Ybs=transform(Ybs,target)
    Ytargetbspred=Clarity_Predict(Ybs,clearned)
    return(Clarity_Extract(Ytargetbspred))
}


###############################
#' @title Compute Residuals or Persistences on bootstrapped data without scaling
#'
#' @description
#' Given a resampled matrix Ybs which was originally used to make a Clarity or ClarityScan object clearned, we:
#' 1. Use the provided clearned object to predict the bootstrapped data
#' 2. Compute the residual structure of this prediction
#'
#' We predict the target with \code{\link{Clarity_Predict}}.
#' 
#' @param clearned either a ClarityScan or Clarity object. If a Clarity object, then residuals will be bootstrapped. If a ClarityScan object then Persistence will be bootstrapped, unless a specific k is given in which case Residuals are again calculated.
#' @param Ybs The bootstrapped matrix for which residuals will be calculated
#' @param summary (default=abs) Summarisation of the bootstrapped value when extracting Residuals. Default: abs, which allows residuals to be easily compared. Use summary=I to return the raw data. See \code{\link{Clarity_Extract}}.
#' @param diag (default=0) diagonal value to be applied when extracting Residuals. See \code{\link{Clarity_Extract}}.
#' @param ... Additional parameters to \code{\link{Clarity_Persistence}} (not required for normal usage)
#' @return A list of residuals if clearned is of class Clarity, or a list of persistences if clearned is of class ClarityScan.
#' @seealso \code{\link{Clarity_Bootstrap}} which is the recommended way to access this code.
#' @export
c_BootstrapValue=function(clearned,Ybs,summary=abs,diag=0,...){
    if(class(clearned)=="ClarityScan"){
        cbs=Clarity_Predict(Ybs,clearned)
        return(Clarity_Persistence(cbs,...))
    }else{
        cbs=Clarity_Predict(Ybs,clearned)
        return(Clarity_Extract(cbs,summary=summary,diag=diag))
    }
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
#' Also note that for normal use you might just want to use \code{\link{plot.Clarity}} or \code{\link{plot.ClarityScan}}, both of which use this function for you.
#'
#' @param testbs A list of bootstrapped matrices as returned from \code{\link{Clarity_Bootstrap}}
#' @param obs Either a Clarity object (for comparing residuals), A ClarityScan object (for comparing persistences) or a matrix of observed values (absolute residuals or persistences, respectively)
#' @param population (default="element") whether the population being compared to is the "element" X[i,j], or "bestcol": the B best values in the "column" X[,i]
#' @param comparison (default=function(x,o){sum(o>x)/(length(x)+1)}) a function to score the observed value relative to the bootstraped values.
#' @param ... Additional parameters to \code{\link{Clarity_Persistence}}, if provided with a ClarityScan object.
#' 
#' @return A matrix of the same shape as obs, containing the result of comparison applied to each element
#' @seealso \code{\link{Clarity_Bootstrap}} for making bootstraps. For plotting, \code{\link{plot.Clarity}} or \code{\link{plot.ClarityScan}} call this function for you.  \code{\link{Clarity_Persistence}} for extracting persistence or  \code{\link{Clarity_Extract}} to extract residuals.
#' @examples
#' \donttest{
#' scan=Clarity_Scan(dataraw) ## Core Clarity
#' predmix=Clarity_Predict(datamix,scan) ## Core prediction
#' 
#' ## Bootstrap persistences:
#' scanbootstrap=Clarity_Bootstrap(scan,target=datamix,D=datarawD)
#' ## Extract observed persistences
#' P=Clarity_Persistence(predmix)
#' ## Compute pvalues
#' pvals=Clarity_BScompare(scanbootstrap,P)
#' ## pvals is a matrix of dimension N by K
#' signif=pvals<0.01
#' ## signif is a logical matrix of dimension N by K
#' 
#' ## Similarly for residuals:
#' ## Bootstrap residuals
#' k10bootstrap=Clarity_Bootstrap(Clarity_Extract(scan,10),
#'                                target=datamix,
#'                                D=datarawD)
#' ## Extract observed residuals
#' k10residuals=Clarity_Extract(predmix,k=10)
#' ## Compute Pvals
#' residualpvals=Clarity_BScompare(k10bootstrap,k10residuals,
#'                         population="bestcol")
#' residualsignif=residualpvals<0.01
#'
#' }
#' @export
Clarity_BScompare=function(testbs,obs,
                           population="element",
                           comparison=function(x,o){sum(c(x,o)>=o)/(length(x)+1)},
                           ... ) {
    if(class(obs)=="Clarity"){
        obs=Clarity_Extract(obs,diag=0)
    }else if(class(obs)=="ClarityScan"){
        obs=Clarity_Persistence(obs,...)
    }
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
            b1=utils::head(sort(b1,decreasing=TRUE),length(testbs))
            sapply(1:dim(tobs)[2],function(j){
                comparison(b1,tobs[i,j])
            })
        }))
    }else(stop("Invalid population in Clarity_BScompare"))
    if(class(ret)=="matrix") {
        rownames(ret)=rownames(obs)
    }else{
        names(ret)=rownames(obs)
    }
    return(ret)
}


###############################
#' @title Perform a Procrustes transformation of one matrix into another
#'
#' @description
#'
#' Procrust the matrix provided in the first argument, into the target provided in the second.  Performs Translation and Dilation.
#' @param Ystart the matrix to transform
#' @param Ytarget The target matrix
#' @return A matrix which is the best fitting transformation of Ystart to Ytarget
#' @export
#' @seealso \code{\link{Clarity_Bootstrap}}, which uses this code.
Clarity_Procrust=function(Ystart,Ytarget){
    if (nrow(Ytarget) != nrow(Ystart)) {
        stop("Dimensions differ between Ystart and Ytarget")
    }
    if (ncol(Ytarget) != ncol(Ystart)) {
        stop("Dimensions differ between Ystart and Ytarget")
    }
    n <- nrow(Ystart)
    m <- ncol(Ystart)
    ## Translation
    J <- diag(n) - 1/n * matrix(1, n, n)
    ## 
    C <- t(Ytarget) %*% J %*% Ystart
    svd.out <- svd(C)
    R <- svd.out$v %*% t(svd.out$u)
    s <- 1
    ## Dilation
    mat1 <- t(Ytarget) %*% J %*% Ystart %*% R
    mat2 <- t(Ystart) %*% J %*% Ystart
    s.numer <- 0
    s.denom <- 0
    for (i in 1:m) {
        s.numer <- s.numer + mat1[i, i]
        s.denom <- s.denom + mat2[i, i]
    }
    s <- s.numer/s.denom
    ## 
    tt <- matrix(0, m, 1)
    ## Translation
    tt <- 1/n * t(Ytarget - s * Ystart %*% R) %*% matrix(1, n, 1)
    ret=s *Ystart %*% R + matrix(tt, n, m, byrow = TRUE)
    ret
}
