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
        if(is.numeric(what)) return(cscan$scan[[what]])
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
Clarity_Persistence=function(clist,f="RowSumsSquared",what="Yresid"){
    if(class(clist)=="ClarityScan") return(Clarity_Persistence(Clarity_Extract(clist,what)))
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
#' @title Diagonal standardization for similarity/dissimilarity matrices
#'
#' @description
#' Complete the diagonal with a sensible approximation that will work for SVD and significance testing. 
#'
#' 
#' @param x An N by N matrix
#' @param method (default="detect") standardization method; the default detects whether a similarity or dissimilarity is provided and uses the max/min respectively. Can be manually set to the function min,max or other operation on a row,column tuple.
#' 
#' @return An N by N matrix of the distances
#' 
#' @export
#' 
c_fixdiagonals=function (x,method="detect"){
    if(dim(x)[1]!=dim(x)[2]) stop("Error: c_fixdiagonals requires a square matrix.")
    if(class(method)!="function"){
        test=sapply(1:dim(x)[1],function(i){
            c(x[i,i],mean(c(x[i,-i],x[-i,i])))
        })
        if(mean(test[1,])<mean(test[2,])){
            method=min
        }else{
            method=max
        }
    }
    ## Distance with the diagonal set to the next-lowest value
    for(i in 1:dim(x)[1]) x[i,i]=method(x[i,-i],x[-i,i])
    x
}

###############################
#' @title Standard Distance function for Clarity
#'
#' @description
#' Computes the euclidean distance between all points using as.matrix(dist(x)), using \code{\link{c_fixdiagonals}} to enforce the diagonal to take the minimum value of the rest of each row.
#'
#' Properly, the diagonal should be be the distance that would be achieved by a replicate of the data x, but this is typically unavailable without making assumptions on the distribution of the features. Using the minimum is typically appropriate; see simulationToDistance in the package ClaritySim for details and the help from that package for a demonstration.
#' 
#' @param x An N by K matrix of N subjects observed at K features
#' 
#' @return An N by N matrix of the distances
#' @seealso \code{\link{Clarity_Compare}}, which takes an argument distfn which can be given this or another function that generates a similarity from feature data.
#' 
#' @export
#' 
c_dist=function (x){
    ## Distance with the diagonal set to the next-lowest value
    d=as.matrix(stats::dist(x))
    return(c_fixdiagonals(d,min))
}

###############################
#' @title Cross validatation method to estimate the number of PCs associated with structure in a Ylist
#'
#' @description
#' From a list of pairs of resampled matrices with the same distribution, use the first to predict the second. The prediction error is then used to form an estimate of the number of dimensions that are increasing predictive power.
#' 
#' This is summarised in 2 ways. Firstly, a point estimate of kmax: the minimum of the average prediction error. Secondly, a soft estimate of Kmax: the proportion of samples that have a minimum prediction error larger than a given K.
#' 
#' @param datalist A list of pairs of N by N matrices, each independently sampled from the same data source as described in \code{\link{Clarity_Predict}}.
#' @param kmax (default=NULL, meaning use the dimension of the data)
#' @param ... additional parameters to \code{link{Clarity_Scan}}.
#' 
#' @return An list containing:
#' \itemize{
#' \item min.score : the point estimate of kmax, obtained from the minimum average prediction error.
#' \item scores : the average prediction error at each complexity \code{K}.
#' \item softk : For each complexity \code{K}, the probability that the estimated \code{K'} is at least as large as \code{K}.
#' }
#' @seealso \code{\link{Clarity_Compare}}, which takes an argument distfn which can be given this or another function that generates a similarity from feature data.
#' 
#' @export
#' 
c_crossvalidate=function(datalist,kmax=NULL,...){
    if(is.null(kmax)){
        d=datalist[[1]]
        kmax=dim(d[[1]])[[1]]-1
    }
    scoremat=sapply(datalist,function(d){
        cp=Clarity_Scan(d[[1]],v=F,kmax=kmax,...)
        scores=sapply(1:length(cp$scan),
                      function(i)sum((cp$scan[[i]]$prediction-d[[2]])^2))
        scores
    })
    avescores=rowMeans(scoremat)
    scoredist=ecdf(apply(scoremat,2,which.min))
    scores=scoredist(1:dim(scoremat)[1])
    return(list(min.score=which.min(avescores), scores=avescores,softk=1-scores))
}


###############################
#' @title Compare a new dissimilarity matrix to another via the learned Clarity model
#'
#' @description
#'
#' Compare a similarity matrix \code{Ynew} to the Clarity model \code{clearned} (learned on \code{Y}) by a) transforming Ynew to maximally look like Y; b) computing persistences or residuals; c) optionally, comparing these to replicated data that produce "alternative versions" of Y.
#'
#' To skip p-value computation, only \code{Ynew} and \code{clearned} are required. In this case, the only difference from Clarity_Predict is that \code{Ynew} is first transformed with \code{transform} to "look like" the original data.
#'
#' To compute p-values, provide additionally a) the original data feature matrix \code{Dnew} and its associated comparison function \code{distfn}, or b) a list of replicated Y (dis)similarity matrices \code{Ylist}. Each replicate must contain a list of three matrices: [Y0,Y1,Y2] where Y0 and Y1 are ferom disjoint samples of the original data, and Y2 is sampled from the new data Dnew similarly to Y1. See \code{\link{c_SplitData}} for how this is done on feature matrices.
#'
#' Some additional checking is done to try to ensure that sane p-values can be computed.
#' 
#' @param clearned A Clarity or ClarityScan object, from which inference will be drawn. If a ClarityScan object and \code{k} is not set, persistences are examined; otherwise residuals are examined.
#' @param Ynew (default=NULL) The new (dis)similarity matrix to be compared. Can instead provide \code{Dnew} and \code{Ynew} will be computed.
#' @param Dnew (default=NULL) Optionally, the new data to be compared. Will be used to construct \code{Ynew} (with \code{distfn}) if that is not provided; otherwise it will be checked for consistency.
#' @param H0 (default="structure") The null hypothesis to be tested. The default test is that the "structure" has changed. Can also be set to "value" in which case no rotation to \code{Y} is done, or "scale" in which case the matrices are rescaled to match by mean value only. (Acts only to set the \code{transform}.)
#' @param D (default=NULL) The *original* data feature matrix. If provided, the features will be resampled to generate replicated datasets and hence (with \code{distfn}) a list of resampled dis/similarity matrices \code{Ylist}.
#' @param Ylist (default=NULL) A list of replicated dis/similarity matrices. Provide them this way if it is more convenient.
#' @param distfn (default=c_dist) Distance function to apply to \code{Dnew} and replicates of \code{D} to produce \code{Ynew} and \code{Ylist}.
#' @param checktol (default=\code{sqrt(.Machine$double.eps)}) tolerance for checking that \code{Ynew} matches that computed from \code{Dnew}.
#' @param transform (default=NULL) Transformation to apply to \code{Ynew} and replicates \code{Ylist} to match them to \code{Y}. This is set by \code{H0} if provided, but can be one of \code{\link{c_Procrust}} (H0="structure") \code{\link{c_MeanScale}} (H0="scale") or \code{\link{c_Identity}} (H0="value")
#' @param comparison (default=\code{\link{c_CompareNormal}}) the comparison between observed and resampled values to compute p-values, controlling the alternative hypothesis H1. Choices are one-tailed tests (\code{\link{c_CompareNormal}}, \code{\link{c_CompareEmpirical}}) looking for large residuals/persistences, or two-tailed tests (\code{\link{c_CompareNormalTwoTailed}} and \code{\link{c_CompareEmpiricalTwoTailed}}) also looking for small values.
#' @param k (default=NULL) If \code{clearned} is a ClarityScan object, the Clarity object with this object is extracted and used instead.
#' @param nbs (defaul=100) Number of resamples to perform, if \code{D} is provided to construct replicate \code{Ylist}.
#' @param verbose (default TRUE) Verbosity level.
#' @param ... Additional arguments to \code{\link{Clarity_Scan}} or \code{\link{Clarity_fixedK}}.
#' 
#' @return If p-values are not computed, an object of the same class as \code{clearned} predicting the new data. Otherwise, a \code{ClaritySignificance} or \code{ClarityScanSignificance} object (if \code{clearned} was a Clarity or ClarityScan object respectively). This is a list consisting of:
#' \itemize{
#' \item Ynew: the new dis/similarity matrix
#' \item Yobs: the transformed new dis/similarity matrix
#' \item pred: the \code{Clarity}/\code{ClarityScan} object predicting \code{Yobs} from \code{clearned}
#' \item splitlist: a list of resampled results as returned from \code{\link{c_SplitComparison}}
#' \item pvals: the computed pvalues
#' \item H0: the provided H0
#' \item transform: the applied transform
#' \item distfn: the applied distfn
#' }
#' @seealso \code{\link{Clarity_Scan}} for computing \code{clearned}; \code{\link{Clarity_Predict}} for prediction; \code{\link{Clarity_Compare}} for assessing uncertainty; \code{\link{plot.Clarity}} and \code{\link{plot.ClarityScan}}.
#' @export
#'
Clarity_Compare <- function(clearned,
                            D=NULL,
                            Dnew=NULL,
                            H0="value",
                            Ylist=NULL,
                            Ynew=NULL,
                            distfn=c_dist,
                            checktol=sqrt(.Machine$double.eps),
                            transform=NULL,
                            comparison=c_CompareNormal,
                            k=NULL,
                            nbs=100,
                            verbose=TRUE,
                            ...
                            ){

### Determine the mode we are running in
    if(is.null(D) && all(is.null(Ylist))){
        mode="predict"
        if(verbose) cat("Using prediction mode as not supplied with Ylist or D to perform resample\n")
    }else{
        mode="resample"
        if(verbose) cat("Using resample mode with supplied data\n")
    }
    
### Determine the type of transform
    if(is.null(transform)){
        if(H0=="structure") {
            transform=c_Procrust
            if(verbose) cat("Using H0=structure, transform=Clarity_Procrust\n")
        }else if(H0=="scale"){
            transform=c_MeanScale
            if(verbose) cat("Using H0=structure, transform=Clarity_Procrust\n")
        }else if(H0=="value"){
            transform=c_Identity
            if(verbose) cat("Using H0=structure, transform=Clarity_Procrust\n")
        }else{
            stop(paste("Error: Invalid H0 :",H0))
        }
    }else{
        if(class(transform)!="function") stop("Error: transform must be a function that accepts two variables, the initial and the target for transformation. See help(\"ClarityCompare\".")
        if(verbose) cat("Using custom H0 with provided transform\n")
    }

### Check data type
    if(!is.null(k) && class(clearned)=="ClarityScan"){
        clearned=clearned$scan[[k]]
        if(verbose) cat(paste0("Extracting Clarity object at requested k=",k," into clearned.\n"))
    }
    
    if(class(clearned)=="ClarityScan"){
        if(verbose) cat("Working with Persistences as provided with a ClarityScan object in clearned.\n")
        kmax=clearned$kmax
        k=NULL
    }else if (class(clearned)=="Clarity"){
        k=clearned$k
        kmax=NULL
        if(verbose) cat(paste0("Working with Residuals as provided with a Clarity object at k=",k,"in clearned.\n"))
    }else{
        stop("Error: Invalid class of clearned object. Must be Clarity or ClarityScan")
    }


### Get Ynew from Dnew if required 
    if(is.null(Ynew)){
        if(is.null(Dnew)) stop("Error: Must provide either Ynew or Dnew")
        if(verbose) cat("Computing Ynew from provided Dnew and distfn.\n")
        Ynew=distfn(Dnew)
    }else if(!is.null(Dnew)){
        testY=distfn(Dnew)
        if(verbose) cat("Using provided Ynew.\n")
        if(!all.equal(Ynew,testY,tolerance=checktol)){
            stop("Error: Provided Ynew and the Y computed from Ynew do not match. This likely means that you have not correctly specified the distfn that you used. If this is a numerical rounding error, increase checktol.")
        }
    }else {
        warning("Could not check that Ynew was created with the same procedure as will be used to perform resampling.")
    }

###    Check the diagonals are structured in Ynew as in Yobs
    if(all(diag(clearned$Y)==0)){
        if(!all(diag(Ynew)==0)){
           warning("Not all Ynew diagonals are zero, whereas the Y from Clarity are?")
        }else{
           warning("The diagonals of your data are zero. In some cases this can prevent correct transformation between datasets resulting is spurious significant structure. Consider revising or using c_fixdiagonals.")
        }
    }else if(!all(diag(clearned$Y)==0)){
        if(all(diag(Ynew)==0)){
           warning("Not all Y from Clarity diagonals are zero, whereas the Ynew are?")
        }
    }

    YnewT=transform(Ynew,clearned$Y)
    pred=Clarity_Predict(YnewT,clearned)
    
    ### Construct the list of resampled data
    if(mode=="resample"){
        if(verbose) cat("Computing resampled data....\n")
        if(all(is.null(Ylist))){
            Ylist=c_SplitData(D,
                              Dnew,
                              nbs=nbs,
                              distfn=distfn)
        }else{
            tlens=sapply(Ylist,length)
            if(any(tlens!=3)) stop("Ylist must be a list, for which each entry is a list of three matrices; Y(resampled from D), Y(resampled from D without reuse of data), and Y(resampled from Dnew)")
        }

### Perform the comparison on each pair of data splits
        if(verbose) cat("Performing resampled predictions....\n")
        csplitlist=c_SplitComparison(Ylist,transform=transform,
                                     k=k,kmax=kmax,verbose=verbose,...)

### Compute p-values
        if(verbose) cat("Performing cross validation to estimate the maximum complexity.\n")
        cv=c_crossvalidate(Ylist)
        if(verbose) cat("Computing pvalues from resampling procedure.\n")
        pvals=c_Scoresplit(csplitlist)
        if(class(pred)=="Clarity"){
            rownames(pvals)=colnames(pvals)=rownames(pred$Y)
        }else if(class(pred)=="ClarityScan"){
            rownames(pvals)=rownames(pred$Y)
        }
        ret=list(
            Ynew=Ynew,
            Yobs=clearned$Y,
            YnewT=YnewT,
            Ylist=Ylist,
            pred=pred,
            splitlist=csplitlist,
            pvals=pvals,
            H0=H0,
            transform=transform,
            distfn=distfn,
            k=k
        )
        if(class(clearned)=="Clarity") class(ret)="ClaritySignificance"
        if(class(clearned)=="ClarityScan") class(ret)="ClarityScanSignificance"
    }else{
        ret=pred
    }
    
    return(ret)
}


###############################
#' @title Emprical estimator for comparing distributions
#'
#' @description
#' the probability that a sample of a distribution (x2) is less than or equal to the sample of another distribution (x1), using their emprical distribution functions. A regularisation is included, to prevent 0s.
#'
#' Probability that x2>=x1 is
#' sum_{i=1}^N pr(x1=x) pr(x2<=x)
#' = (1/(reg+length(x))) e2(c(reg,x1))
#'
#' @param x1 A set of samples from the reference distribution
#' @param x2 A set of samples from the query distribution
#' @return the probability that a sample from the query is larger than the reference.
#' @export
#' @seealso \code{\link{Clarity_Compare}}, the recommended way to access this code.
#' @examples
#' c_px2_lt_x1(0,1) # 0.5 due to regularization
#' c_px2_lt_x1(0,1,0) # 0
#' c_px2_lt_x1(0:2,1:3,0) # 1/3
c_px2_lt_x1<-function(x1,x2,reg=1){
    d1=sort(x1)
    d2=sort(x2)
    e2=stats::ecdf(d2)
    ##  sum(e2(d1)/(1+length(d1)))
    sum(c(reg,e2(d1))/(reg+length(d1)))
}

###############################
#' @title Comparison for lists of pairs of matrices
#'
#' @description
#' Perform an element-wise comparison on the set of elements at a given location in a matrix, computing:
#' \code{comparison({X[i,j]}_{1..B},{Y[i,j]}_{1..B})}
#' i.e comparing the set of [i,j]-th elements of Y to the same set for X.
#'
#' The intended use case is when X and Y are each generated by some scoring function for resamples of a reference dataset (X) and a query dataset (Y). The default uses \code{\link{c_px2_gt_x1}} to calculate the probability that Y[i,j] is larger than X[i,j].
#' 
#' @param scorelist A list in which each element is a length 2 list, consisting of [X,Y] where X and Y are each matrices of equal dimension.
#' @param comparison (default=c_px2_gt_x1) comparison function to evaluate.  \code{function(d1,d2)=wilcox.test(d1,d2)$p.value} and similar distributional tests are more stringent but are not appropriate for resampled data.
#' 
#' @return the probability that a sample from the query is larger than the reference.
#' @export
c_Scoresplit=function(scorelist,comparison=c_px2_lt_x1){
    size1=dim(scorelist[[1]][[1]])[1]
    size2=dim(scorelist[[1]][[1]])[2]
    ret=t(sapply(1:size1,function(i){
        sapply(1:size2,function(j){
            d1=sapply(1:length(scorelist),function(rep) scorelist[[rep]][[1]][i,j])
            d2=sapply(1:length(scorelist),function(rep) scorelist[[rep]][[2]][i,j])
            comparison(d1,d2)
        })
    }))
    ret
}


###############################
#' @title Split two datasets into three parts by downsampling
#'
#' @description
#' Downsample columns of D1 and D2 into three components:
#'   D1'': Random half of the columns (rounded up) of D1
#'   D1' : The other half of the columns (rounded down) of D1
#'   D2' : Random half of the columns (rounded down) of D2
#' distfn is evaluated on each to return (dis)similarities.
#'
#' @param D1 reference matrix consisting N items with L1 features 
#' @param D2 reference matrix consisting N items with L2 features
#' @param distfn (default: \code{\link{c_dist}}) summary to apply, which should return a matrix of size N times N
#' @param nbs (default: 100) Number of times to repeat this process (set to NULL to do it only once and return a single resample list  [Y1'',Y1',Y2']
#' @return A list of length nbs, each of which is a list containing [Y1'',Y1',Y2']
#' @export
#' @seealso \code{\link{Clarity_Compare}}, the recommended way to access this code.
c_SplitData=function(D1,D2,distfn=c_dist,nbs=100){
    if(is.null(nbs)){
        xbs1 = sample(1:dim(D1)[2], floor(dim(D1)[2]/2), replace = FALSE)
        xbs2 = sample(1:dim(D2)[2], floor(dim(D2)[2]/2), replace = FALSE)
        D1ref = D1[, -xbs1, drop = FALSE]
        D1sample = D1[, xbs1, drop = FALSE]
        D2sample = D2[, xbs2, drop = FALSE]
        Ybs = list(Y1=distfn(D1ref),
                   Y2=distfn(D1sample),
                   Y3=distfn(D2sample))
        return(Ybs)
    }else{
        return(lapply(1:nbs,function(i)c_SplitData(D1,D2,distfn=distfn,nbs=NULL)))
    }
}

###############################
#' @title Perform a comparison on a list of resampled matrices
#'
#' @description
#' Given a list of [Y1'', Y1', Y2'] matrix triples, learn a Clarity model on Y1'' , and predict Y1' and Y2' from that model. Do this for each triple and report the results as a pair of distributions, that for Y1' and for Y2'.
#'
#' @param Ydata1 Either a list of matrix triples, or Y1'', the matrix to use to construct the Clarity model
#' @param Ydata2 (default=NULL): If NULL, we assume Ydata1 is a list of matrix triples. If non-null, the matrix Y1', i.e. replicated data from the same distribution as Ydata1
#' @param Yalt (default=NULL): Should be NULL if Ydata2 is NULL. If non-null, should be Y2', that is a sampled alternative (dis)similarity matrix to be queried.
#' @param transform (default=\code{link{c_Procrust}}) a function taking two values, the first of which is transformed into the second. Use =\code{link{c_Identity}} to disable.
#' @param k (default=NULL) If null, perform \code{\link{Clarity_Scan}} prediction for Persistences. Otherwise perform \code{\link{Clarity_fixedK}} at the specified k.
#' @param verbose (default=FALSE) verbosity if operating over a list of matrix triples. 
#' @param ... Additional parameters to Clarity. The kmax parameter for  \code{\link{Clarity_Scan}} could be important.
#' @return If provided with a list of matrix triples, return a list of matrix pairs giving the summary of the Clarity model predicting [Ydata2,Yalt]. If provided with three matrices, return two matrix predicting [Ydata2,Yalt]. For fixed k, this is a matrix of Residuals from \code{\link{Clarity_Extract}}; otherwise a matrix of Persistences from \code{\link{Clarity_Persistence}}.
#' @export
#' @seealso \code{\link{Clarity_Compare}}, the recommended way to access this code.
c_SplitComparison=function(Ydata1,
                           Ydata2=NULL,
                           Yalt=NULL,
                           transform=c_Procrust,
                           k=NULL,
                           verbose=FALSE,
                           ...){
    if(is.null(Ydata2)){## We will iterate over a list of matrix triples
        Ylist=Ydata1
        nbs=length(Ylist)
        bspers = lapply(1:nbs, function(rep) {
            if (verbose) 
                cat(paste("... replication", rep, "of", nbs, "\n"))
            
            return(c_SplitComparison(Ylist[[rep]][[1]],
                                     Ylist[[rep]][[2]],
                                     Ylist[[rep]][[3]],
                                     transform=transform,
                                     k=k,...))
        })
        return(bspers)
    }else{ ## We have been given a single matrix triple to compare
        Ydata2=transform(Ydata2,Ydata1)
        Yalt=transform(Yalt,Ydata1)
        if(is.null(k)){ ## Use ClarityScan
            clearned=Clarity_Scan(Ydata1,verbose=FALSE,...)
            cdata=Clarity_Predict(Ydata2,clearned)
            cpred=Clarity_Predict(Yalt,clearned)
            
            pdata=Clarity_Persistence(cdata)
            ppred=Clarity_Persistence(cpred)
        }else{ ## Use Clarity at the specified k
            clearned=Clarity_fixedK(Ydata1,k,verbose=FALSE,...)
            cdata=Clarity_Predict(Ydata2,clearned)
            cpred=Clarity_Predict(Yalt,clearned)
            
            pdata=Clarity_Extract(cdata)
            ppred=Clarity_Extract(cpred)
        }
        return(list(pdata=pdata,ppred=ppred))
    }
}
