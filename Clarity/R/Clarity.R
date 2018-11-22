
###############################
#' @title Initialise A and X in Clarity using Kmeans
#' @description
#' This function takes a non-negative similarity matrix Y of dimension N by N, 
#' and a desired number of clusters K,
#' and uses K-means to estimate them.
#' 
#' It then "noises" these clusters by retaining a membership fraction a of the chosen cluster and distributing membership (1-a) amongst the other K-1 clusters, with additional multiplicative noise N(1,sigma).
#'
#' @param Y An N by N non-negative matrix of similarities
#' @param K number of clusters, passed to \code{\link{kmeans}}
#' @param a (default 0.9) Fraction of mixture to be retained in the chosen cluster
#' @param sigma (default 0.0001) Multiplicative noise for A
#' @param sigmaY (default 1e-5) Multiplicative noise for Y, to prevent co-linearity
#' @param nstart (default 1000) Number of restarts, passed to \code{\link{kmeans}}
#' @param iter.max (default 50) Number of iterations, passed to \code{\link{kmeans}}
#'  
#' @keywords mixture
#' @return A list containing the following:
#' \itemize{
#' \item Y An N by N non-negative matrix of similarities
#' \item A An N by K non-negative matrix of mixture proportions
#' \item X A K by K non-negative matrix of latent weights
#' }
#' @seealso \code{\link{Clarity_Scan}} is the recommended interface for performing  inference on a representation of a single similarity matrix with Clarity.
#' @export
#' @examples

c_initKmeans=function(Y,K,a=0.9,sigma=0.0001,sigmaY=1e-5,nstart=1000,iter.max=50){
    if(dim(Y)[1]==K){
        A=diag(K)
        H=Y
        X=c_updateX(A,H,Y)
        return(list(Y=Y,A=A,X=X))
    }
    if(K==1){
        A=matrix(1,nrow=dim(Y)[1],ncol=1)
        H=1
        X=c_updateX(A,H,Y)
        return(list(Y=Y,A=A,X=X))
    }
    V=Y+stats::runif(dim(Y)[1]*dim(Y)[2],0,mean(abs(Y))*sigmaY)
    test=stats::kmeans(Y,K,nstart=nstart,iter.max=iter.max)
    H=test$centers
    A=t(apply(V,1,function(x){
        tdists=colSums((x-t(H))^2)
        ret=rep((1-a)/(dim(H)[1]-1),dim(H)[1])
        ret[which.min(tdists)]=a
        ret=ret*stats::rnorm(length(ret),1,sigma)
        ret/sum(ret)
    }))
    X=c_updateX(A,H,Y)
    return(list(Y=Y,A=A,X=X))
}


###############################
#' @title Multplicative update for A in Clarity
#'
#' @description
#' This function improves the fit of
#' 
#'      Y = || A X A^T ||_2
#' 
#' by updating A. X and Y are assumed to be positive, and A is positive with rows summing to 1. A has dimension N by K.
#'
#' Seeking to solve
#' 
#'             Y = A X t(A)
#' 
#' we minimise F = (1/2) || Y - A X t(A) ||_2
#' 
#' The calculation leads to
#' 
#' A_{t+1} = A [ t(Y) A X + Y A t(X) ]
#' 
#'             ---------------------------------
#' 
#'  [ A X t(A) A t(X) + A t(X) t(A) A X ]
#'
#'
#' @param A An N by K non-negative matrix of mixture proportions
#' @param X A K by K non-negative matrix of latent weights
#' @param Y An N by N non-negative matrix of similarities
#' @param fixed IGNORE
#' @param norm Whether to return the updated A normalised to have row sum 1. Default: TRUE.
#'
#' @keywords mixture
#' @return A non-negative matrix of dimension N by K
#' @export
#' @examples
c_updateA<-function(A,X,Y,fixed=NULL,norm=TRUE) {
    ## This is probably the most efficient way to implement this, based on the choices of ordering available
    tAA=(t(A) %*% A) 
    num = t(Y) %*% A %*% X + Y %*% A %*% t(X)
    denom = A %*% X %*% tAA %*% t(X) + A %*% t(X) %*% tAA %*% X
    ret=A * (num / denom)
    if(norm)ret=ret/rowSums(ret)
    ret
}

###############################
#' @title Exact update for X in Clarity
#'
#' @description
#' This function improves the fit of
#'      Y = || A X A^T ||_2
#' by updating X. X and Y are assumed to be positive, and A is positive with rows summing to 1. A has dimension N by K.
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
#' @examples
c_updateX_solve<-function(A,Y){
    tAAinv=try(
        solve(t(A)%*% A),silent=TRUE)
    if(class(tAAinv)=="try-error") return(tAAinv)
    Xsol=tAAinv %*% t(A) %*% Y %*% A %*% tAAinv
    Xsol
}

###############################
#' @title Safe update for X in Clarity
#'
#' @description
#' This function improves the fit of
#'      Y = || A X A^T ||_2
#' by updating X. X and Y are assumed to be positive, and A is positive with rows summing to 1. A has dimension N by K.
#'
#' We use \code{\link{c_updateX_solve}} where possible. However, sometimes X is of low rank and the calculation fails; we fall back on a mulitiplicative update in this case.
#'
#' Seeking to solve 
#'             Y = A X t(A)
#' we minimise F = (1/2) || Y - A X t(A) ||_2 
#' The calculation leads to
#' X_{t+1} = X t(A) Y A / [ t(A) A X t(A) A ]
#'
#' @param A An N by K non-negative matrix of mixture proportions
#' @param X A K by K non-negative matrix of latent weights
#' @param Y An N by N non-negative matrix of similarities
#' 
#' @keywords mixture
#' @return A non-negative matrix of dimension K by K
#' @seealso \code{\link{c_updateA}} for updates to A
#' @export
#' @examples
c_updateX<-function(A,X,Y){
    Xsol=try(c_updateX_solve(A,Y),silent=TRUE)
    if(class(Xsol)=="try-error"){
        tAA=(t(A) %*% A) 
        num = t(A) %*% (Y) %*% A
        denom = tAA %*% X %*% tAA
        Xsol=X * (num / denom)
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
#' @examples
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
#' @examples
c_objectivefunction=function(A,X,Y){
    Yhat=A %*% X %*% t(A) 
    sum((Y - Yhat)^2)
}


###############################
#' @title Fit a Similarity matrix using Clarity
#'
#' @description
#' This function fits the best solution of A and X to
#'      || Y - A X A^T ||_2
#' for a FIXED dimension of the latent space K. Here,
#' A is n N by K non-negative matrix of mixture proportions
#' and
#' X is a K by K non-negative matrix of latent weights
#' 
#' @param Y An N by N non-negative matrix of similarities
#' @param K an integer giving the dimension of the fit to be used
#' @param tmax (default 1e4) how many iterations to run for. maximimum
#' @param matrixdistmin (default 1e-8) stopping criterion for the summed square difference of A and X
#' @param objectivedistmin (default 1e-3) stopping criterion for the increase in objective function
#' @param printskip (default 100) how many iterations between progress updates
#' @param verbose (default TRUE) whether to output any information regarding progress
#' @param updateA (default TRUE) whether to update A; shouldn't use this option
#' @param updateX (default TRUE) whether to update X; shouldn't use this option
#' @param A (default NULL) initial A. Must pass A and X together to specify
#' @param X (default NULL) initial X. Must pass A and X together to specify
#' @param a (default 0.9) a parameter for \code{\link{c_initKmeans}}
#' @param sigma (defailt 0.0001) a parameter for \code{\link{c_initKmeans}}
#' @param sigmaY (defailt 1e-5) a parameter for \code{\link{c_initKmeans}}
#' @param ensurepositive (default TRUE) Whether we restrict to non-negative Y. Note that the algorithm can work with some negative Y as long as it is postitive semi-definite, but this can still lead to convergence problems.
#' @param fixed (efault NULL) which columns to keep structurally fixed. Not yet used.
#' 
#' @keywords mixture
#' @return A list of class "Clarity" containing the following:
#' \itemize{
#' \item Y, the provided non-negative N by N similarity matrix
#' \item A, the inferred non-negative N by K mixture matrix
#' \item X, the inferred non-negative K by K latent matrix
#' \item calX, the inferred non-negative K by K 'naturally scaled' latent matrix
#' \item prediction, A X A^T
#' \item C, the matrix C above relating X to calX
#' \item tCinv, the inverse of C^T, used in relating X to calX
#' \item Afirst, the initial estimate of A
#' \item Xfirst, the initial estimate of X
#' \item calXfirst, the initial estimate of calX
#' \item numiters, the number of actually used iterations
#' \item matrixdist, the distance between matrices, a vector of length numiters
#' \item objective, the objective function at the final (A,X;Y)
#' \item objectivedist, the objective function, a vector of length numiters
#' \item matrixdistmin, the matrix change threshold provided
#' \item objectivedistmin, the objective change threshold provided
#' }
#' 
#' @seealso \code{\link{Clarity_Scan}} is the recommended interface for performing  inference on a representation of a single similarity matrix with Clarity.
#' @export
#' @examples
Clarity_fixedK<-function(Y,
                        K=NULL,
                  tmax=1e4, # Maximum number of iterations
                  matrixdistmin=1e-8, # Stopping criterion 
                  objectivedistmin=1e-3, # Stopping criterion 
                  printskip=100, # Report progress every Xth iteration
                  verbose=TRUE, # Do we print to the screen
                  updateA=TRUE, # Do we update A?
                  updateX=TRUE, # Do we update X?
                  A=NULL,
                  X=NULL,
                  a=0.9, # Passed to c_initKmeans 
                  sigma=0.0001, # Passed to c_initKmeans
                  sigmaY=1e-5, # Passed to c_initKmeans
                  ensurepositive=TRUE, # Check for whether Y is strictly non-negative
                  fixed=NULL # Do we fix any columns of A?
                  ) {
    if(ensurepositive){
        if(any(Y<0)) {
            stop("ERROR: Some Y values are negative! You can either a) construct a positive Y, e.g. as.matrix(dist(Y)) or using a different construction, or b) set ensurepositive=FALSE in your call to Clarity. However, this may converge less well.")
        }
    }
    ## Initialisation
    if(is.null(A) || is.null(X)) {
        init=c_initKmeans(Y,K,a=a,sigma=sigma,sigmaY=sigmaY)
        Y = init$Y
        A = init$A
        X = init$X
    }
    ## Record keeping
    matrixdistlist=objectivelist=rep(0,tmax)
    Afirst=A
    Xfirst=X
    calXlist=c_calXfromX(A,X)
    calXfirst=calXlist$calX
    distdiff=dist
    matrixdistlist[1]=sum(A^2) + sum(X^2)
    objectivelist[1]=c_objectivefunction(A,X,Y)
    ## Main loop
    for(t in 2:tmax){
        Alast=A
        Xlast=X
        ## Do the updates
        if(updateA) A=c_updateA(A,X,Y,fixed=fixed)
        if(updateX) X=c_updateX(A,X,Y)

        ## Evaluation metrics
        diff=sum((A-Alast)^2) + sum((X-Xlast)^2) # L2 norm of the change in A and X
        dist=c_objectivefunction(A,X,Y) # Total distance from the objective
        distdiff=objectivelist[t-1]-dist # change in the objective function
        matrixdistlist[t]=diff
        objectivelist[t]=dist
        ## Should we report?
        if((t%%printskip)==0){
            tmpdists=c(last=c_objectivefunction(Alast,Xlast,Y),
                       intermediate=c_objectivefunction(A,Xlast,Y),
                       final=c_objectivefunction(A,X,Y))
            if(verbose) print(paste("Iteration",t,"Difference",format(diff,signif=2),"Objective functions",paste(format(tmpdists,signif=2),collapse=","),"Improvement",format(distdiff,signif=2)))
        }
        ## Should we return to the previous estimate and stop?
        if(objectivelist[t]>objectivelist[t-1]) {
            objectivelist=objectivelist[1:(t-1)]
            matrixdistlist=matrixdistlist[1:(t-1)]
            A=Alast
            X=Xlast
            if(verbose) print(paste("Stopping at iteration",t-1,"due to worsening objective function",dist,", improvement",distdiff,"difference",diff))
            break;
        }
        ## Should we just stop?
        if(diff<matrixdistmin | distdiff<objectivedistmin){
            objectivelist=objectivelist[1:t]
            matrixdistlist=matrixdistlist[1:t]
            if(verbose) print(paste("Stopping at iteration",t,"with objective function",dist,", improvement",distdiff,"difference",diff))
            break;
        }
    }
    ## Final objects
    calXlist=c_calXfromX(A,X)
    prediction=A %*% X %*% t(A) 
    ret=list(Y=Y,A=A,X=X,calX=calXlist$calX,
                prediction=prediction,
                C=calXlist$C,
                tCinv=calXlist$tCinv,
                Afirst=Afirst,
                Xfirst=Xfirst,
                calXfirst=calXfirst,
                numiters=t,
             matrixdist=matrixdistlist,
             objective=utils::tail(objectivelist,1),
                objectivedist=objectivelist,
                matrixdistmin=matrixdistmin,
             objectivedistmin=objectivedistmin)
    class(ret)="Clarity"
    return(ret)
}

###############################
###############################
###############################
###############################

###############################
#' @title Scan over K to learn a list of Clarity objects
#'
#' @description
#' Runs Clarity_fixedK for all K up to some threshold. Returns a list of results which is improved upon by many other functions.
#' 
#' @param Y An N by N non-negative matrix of similarities
#' @param kmax an integer giving the maximum dimension of the fit to be used (2..kmax are used)
#' @param tmax (default 1e4) how many iterations to run for. maximimum
#' @param matrixdistmin (default 1e-8) stopping criterion for the summed square difference of A and X
#' @param objectivedistmin (default 1e-3) stopping criterion for the increase in objective function
#' @param printskip (default 100) how many iterations between progress updates
#' @param verbose (default TRUE) whether to output any information regarding progress
#' @param ... futher parameters to be passed to \code{\link{Clarity_fixedK}}
#' 
#' 
#' @keywords mixture
#' @return A object of class "ClarityScan" as returned by \code{\link{Clarity_Scan}}
#' 
#' @seealso \code{\link{Clarity_Scan}} is the recommended interface for performing  inference on a representation of a single similarity matrix with Clarity.
#' @export
#' @examples
c_ScanOnly=function(Y,kmax=20,tmax=1e4,
                     matrixdistmin=1e-8,objectivedistmin=1e-3,
                    printskip=1000,verbose=TRUE,...){
    klist=1:kmax
    clist=lapply(klist,function(k){
        print(k)
        Clarity_fixedK(Y,K=k,
                      matrixdistmin=matrixdistmin,
                      objectivedistmin=objectivedistmin,
                      printskip=printskip,
                      tmax=tmax,...)
    })
    ret=list(scan=clist,
             objectives=sapply(clist,function(x)x$objective),
             klist=klist,
             Y=Y,
             kmax=kmax)
    class(ret)<-"ClarityScan"
    return(ret)
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
#' @examples
c_listscore=function(clist){
    if(class(clist)!="ClarityScan") stop("clist must be of class ClarityScan as returned by Clarity_Scan")
    sapply(clist$scan,function(x)c_objectivefunction(x$A,x$X,x$Y))
}

###############################
#' @title Make a larger Clarity object by randomly peturbing a smaller one
#'
#' @description
#' Make a Clarity object of size K+1 from one of size K by mixing in a random matrix with the old one
#' 
#' @param Y Non-negative Similiarty matrix being fitted (N by N)
#' @param AatKm1 Non-negative matrix of mixtures (N by K)
#' @param XatKm1 Non-negative matrix of latent structures (K by K)
#' @param alpha (default 0.9) Mixing weight between original matrix and new matrix
#' @param dirichletbeta (default 1) Dirichlet beta for the new random A being mixed with the original
#' 
#' @keywords mixture
#' @return A list containing the following:
#' \itemize{
#' \item X A K by K non-negative matrix of latent weights
#' \item A An N by K non-negative matrix of mixture proportions
#' \item dist The objective function evaluated at the new A and X
#' }
#' 
#' @export
#' 
c_IncreaseK<-function(Y,AatKm1,XatKm1,alpha=0.9,dirichletbeta=1){
    ## Pick an initial point with increased dimension as a random pertubation from a run at K-1
    km1=dim(AatKm1)[2]
    k=km1+1
    n=dim(AatKm1)[1]
    trandom=rdirichlet(n,rep(dirichletbeta,length=k))
    testA=alpha * cbind(AatKm1,0) + (1-alpha)*trandom
    X0=cbind(rbind(XatKm1,0),0)
    Xrandom=rdirichlet(k,rep(1,length=k))
    testX=alpha*X0 + (1-alpha)*Xrandom # In case we have to fall back to sequential update
    testX=c_updateX(testA,testX,Y) 
    list(X=testX,A=testA,dist=c_objectivefunction(testA,testX,Y))
}

###############################
#' @title Make a smaller Clarity object by randomly peturbing a larger one
#'
#' @description
#' Make a Clarity object of size K-1 from one of size K-1 by removing a random column from of size K
#' 
#' @param Y Non-negative Similiarty matrix being fitted (N by N)
#' @param A Non-negative mixture matrix 
#' @param X Non-negative latent similarity matrix
#' 
#' @keywords mixture
#' @return A list containing the following:
#' \itemize{
#' \item X A K by K non-negative matrix of latent weights
#' \item A An N by K non-negative matrix of mixture proportions
#' \item dist The objective function evaluated at the new A and X
#' }
#' 
#' @export
#' 
c_DecreaseK<-function(Y,A,X){
    ## Pick an initial point with reduced dimension as a random dropping of a run at K+1
    if(dim(A)[2]<=1) stop("Cannot reduce matrix size!")
    tdrop=sample(1:dim(A)[2],1)
    testA=A[,-tdrop,drop=FALSE]
    testA=testA/rowSums(testA)
    if(any(is.na(testA))) testA[is.na(testA)[,1],]=1/dim(testA)[2]
    testX=X[-tdrop,-tdrop,drop=FALSE] # In case we can't solve for X
    testX=c_updateX(testA,testX,Y) 
    list(X=testX,A=testA,dist=c_objectivefunction(testA,testX,Y))
}


###############################
#' @title Update a Clarity estimate by scanning over increasing K
#'
#' @description
#' Takes a list of Clarity results at different K and update them all, using estimates at lower K to improve those at higher K.
#' 
#' @param clist A list of Clarity objects as returned by \code{\link{Clarity_Scan}}
#' @param verbose Whether to output information about progress. Note that Clarity_fixedK is always run with verbose=FALSE
#' @param tmax tmax as passed to \code{\link{Clarity_fixedK}}
#' @param matrixdistmin (default NULL meaning derive from clist) Passed to \code{\link{Clarity_fixedK}}.
#' @param objectivedistmin (default NULL meaning derive from clist) Passed to \code{\link{Clarity_fixedK}}.
#' @param ... futher parameters to be passed to \code{\link{Clarity_fixedK}}
#' 
#' @keywords mixture
#' @return A object of class "ClarityScan" as returned by \code{\link{Clarity_Scan}}
#' 
#' @seealso \code{\link{Clarity_Scan}} is the recommended interface for performing  inference on a representation of a single similarity matrix with Clarity.
#' @export
#' @examples
c_ForwardStep<-function(clist,verbose=TRUE,tmax=10000,matrixdistmin=NULL,objectivedistmin=NULL,...){
    ## Update states at K+1 using the state at K, and choose the best one to retain
    ## Report the new updated list
    if(class(clist)!="ClarityScan") stop("clist must be of class ClarityScan as returned by Clarity_Scan")
    scan=clist$scan
    if(is.null(matrixdistmin))matrixdistmin=scan[[1]]$matrixdistmin
    if(is.null(objectivedistmin))objectivedistmin=scan[[1]]$objectivedistmin

    for(i in 1:(length(scan)-1)){
        ## Continue a previous run
        ret0=Clarity_fixedK(scan[[i+1]]$Y,A=scan[[i+1]]$A,X=scan[[i+1]]$X,
                     tmax=tmax,
                     matrixdistmin=matrixdistmin,
                     objectivedistmin=objectivedistmin,
                     verbose=F,...)
        ## Do a new run at an increased K
        init=c_IncreaseK(scan[[i]]$Y,scan[[i]]$A,
                                scan[[i]]$X,alpha=0.9,dirichletbeta=1)
        ret=Clarity_fixedK(scan[[i]]$Y,A=init$A,X=init$X,tmax=tmax,
                     matrixdistmin=matrixdistmin,
                     objectivedistmin=objectivedistmin,
                     verbose=F,...)
        scores=c(orig=scan[[i+1]]$objective,
                 base=ret0$objective,
                 seq=ret$objective
                 )
        if(scores[1]==min(scores)){
            if(verbose) print(paste("list entry",i+1,"rejected changes"))
        }else if(scores[3]==min(scores)){ ## New run is better
            if(verbose) print(paste("list entry",i+1,"Using new run with score",
                                    scores[3],"over",scores[1]))
            scan[[i+1]]=ret
        }else if(scores[2]==min(scores)){
            if(verbose) print(paste("list entry",i+1,"Using extended run with score",
                                    scores[2],"over",scores[1]))
            scan[[i+1]]=ret0
        }
    }
    clist$scan=scan
    clist$objectives=sapply(scan,function(x)x$objective)
    clist
}
    
###############################
#' @title Update a Clarity estimate by scanning over decreasing K
#'
#' @description
#' Takes a list of Clarity results at different K and update them all, using estimates at higher K to improve those at lower K.
#' 
#' @param clist A list of Clarity objects as returned by \code{\link{Clarity_Scan}}
#' @param verbose Whether to output information about progress. Note that Clarity_fixedK is always run with verbose=FALSE
#' @param tmax tmax as passed to \code{\link{Clarity_fixedK}}
#' @param matrixdistmin (default NULL meaning derive from clist) Passed to \code{\link{Clarity_fixedK}}.
#' @param objectivedistmin (default NULL meaning derive from clist) Passed to \code{\link{Clarity_fixedK}}.
#' @param ... futher parameters to be passed to \code{\link{Clarity_fixedK}}
#' 
#' @keywords mixture
#' @return A object of class "ClarityScan" as returned by \code{\link{Clarity_Scan}}
#' 
#' @seealso \code{\link{Clarity_Scan}} is the recommended interface for performing  inference on a representation of a single similarity matrix with Clarity.
#' @export
#' @examples
c_BackwardStep<-function(clist,verbose=TRUE,tmax=10000,matrixdistmin=NULL,objectivedistmin=NULL,...){
    ## Update states at K-1 using the state at K, and choose the best one to retain
    ## Report the new updated list
    if(class(clist)!="ClarityScan") stop("clist must be of class ClarityScan as returned by Clarity_Scan")
    if(is.null(matrixdistmin))matrixdistmin=clist$scan[[1]]$matrixdistmin
    if(is.null(objectivedistmin))objectivedistmin=clist$scan[[1]]$objectivedistmin
    scan=clist$scan

    for(i in (length(scan)):2){
        ret0=Clarity_fixedK(scan[[i-1]]$Y,A=scan[[i-1]]$A,X=scan[[i-1]]$X,
                     tmax=tmax,
                     matrixdistmin=matrixdistmin,
                     objectivedistmin=objectivedistmin,
                     verbose=F,...)
        init=c_DecreaseK(scan[[i]]$Y,scan[[i]]$A,
                                         scan[[i]]$X)
        ret=Clarity_fixedK(scan[[i]]$Y,A=init$A,X=init$X,tmax=tmax,
                     matrixdistmin=matrixdistmin,
                     objectivedistmin=objectivedistmin,
                     verbose=F,...)
        scores=c(orig=scan[[i-1]]$objective,
                 base=ret0$objective,
                 seq=ret$objective
                 )
        if(scores[1]==min(scores)){ ## New run is better
            if(verbose) print(paste("Backwards: list entry",i-1,"rejected changes"))
            scan[[i-1]]=ret0
        }else if(scores[3]==min(scores)){ ## New run is better
            if(verbose) print(paste("Backwards: list entry",i-1,"Using new run with score",
                                    scores[3],"over",scores[1]))
            scan[[i-1]]=ret
        }else if(scores[2]==min(scores)){
            if(verbose) print(paste("Backwards: list entry",i-1,"Using extended run with score",
                                    scores[2],"over",scores[1]))
            scan[[i-1]]=ret0
        }
    }
    clist$scan=scan
    clist$objectives=sapply(scan,function(x)x$objective)
    clist
}

###############################
#' @title Update a Clarity estimate by scanning over K 
#'
#' @description
#' Takes a list of Clarity results at different K and update them all, by iteratively trying to improve lower K estimates from higher K, and vice-versa
#' 
#' @param clist A list of Clarity objects as returned by \code{\link{Clarity_Scan}}
#' @param niter (default 10) maximum number of forward-backward iterations
#' @param thresh (default -1) required improvement in objective function across all K before we give up. If negative, then we always try niter iterations. Setting a small, positive value can speen computation in large problems.
#' @param verbose Whether to output information about progress. Note that Clarity_fixedK is always run with verbose=FALSE
#' @param tmax tmax as passed to \code{\link{Clarity_fixedK}}
#' @param matrixdistmin (default NULL meaning derive from clist) Passed to \code{\link{Clarity_fixedK}}.
#' @param objectivedistmin (default NULL meaning derive from clist) Passed to \code{\link{Clarity_fixedK}}.
#' @param ... futher parameters to be passed to \code{\link{Clarity_fixedK}}
#' 
#' @keywords mixture
#' @return A object of class "ClarityScan" as returned by \code{\link{Clarity_Scan}}
#' 
#' @seealso \code{\link{Clarity_Scan}} is the recommended interface for performing  inference on a representation of a single similarity matrix with Clarity.
#' @export
#' @examples
c_ForwardBackward<-function(clist,niter=10,thresh=-1,tmax=10000,verbose=TRUE,matrixdistmin=1e-8,objectivedistmin=0.001,...){
    ## Runs backwards and forwards through the results using neighbouring runs to improve one-another. Does this niter times in both directions. Returns the optimised clist
    if(class(clist)!="ClarityScan") stop("clist must be of class ClarityScan as returned by Clarity_Scan")
    converged=FALSE
    for(i in 1:niter){
        print(paste("Iteration",i))
        score0=c_listscore(clist)
        clist1A=c_ForwardStep(clist,
                                verbose=verbose,tmax=tmax,
                                matrixdistmin=matrixdistmin,
                                objectivedistmin=objectivedistmin,
                                ...)
        score1A=c_listscore(clist1A)
        clist1B=c_BackwardStep(clist1A,
                                 verbose=verbose,tmax=tmax,
                                 matrixdistmin=matrixdistmin,
                                 objectivedistmin=objectivedistmin,
                                 ...)
        score1B=c_listscore(clist1B)
        clist=clist1B
         
        scores=rbind(score0,score1A,score1B)
        if(all(diff(score1B)<=0)
           && (min(rowSums(scores))+thresh>sum(scores[1,])) ) {
            if(verbose)print(paste("Terminating at iteration",i,"due to convergence criteria being met"))
            converged=TRUE
            break;
        }
    }
    if(!converged){
        if(verbose)print(paste("Terminating without clear convergence after trying",niter,"iterations"))
    }
    clist
}

###############################
#' @title Run The full Clarity algorithm
#'
#' @description
#' Learn a mixture model representation of a provided similarity matrix Y using Clarity. You may either choose the maximum number of components to use, kmax, or provide a pre-existing list of 
#' 
#' @param Y An N by N non-negative matrix of similarities
#' @param kmax (default 20) an integer giving the maximum dimension of the fit to be used (2..kmax are used). Ignored if clist is provided.
#' @param clist (default NULL) A list of Clarity objects as returned by \code{\link{Clarity_Scan}}. If not provided, \code{\link{c_ScanOnly}} is used to generate an initial run with the specified kmax.
#' @param niter (default 10) maximum number of forward-backward iterations
#' @param thresh (default 0.1) required improvement in objective function across all K before we give up
#' @param verbose Whether to output information about progress. Note that Clarity_fixedK is always run with verbose=FALSE
#' @param tmax tmax as passed to \code{\link{Clarity_fixedK}}
#' @param matrixdistmin (default NULL meaning derive from clist) Passed to \code{\link{Clarity_fixedK}}.
#' @param objectivedistmin (default NULL meaning derive from clist) Passed to \code{\link{Clarity_fixedK}}.
#' @param ... futher parameters to be passed to \code{\link{Clarity_fixedK}}
#' 
#' @keywords mixture
#' @return An object of class "ClarityScan", which is a list containing the following:
#' \itemize{
#' \item scan A list of length (kmax)-1, each containing an object of class "Clarity" as described in \code{\link{Clarity_fixedK}} (for K=2..kmax)
#' \item objectives A vector of length (kmax)-1, the objective function at each K
#' \item klist The list of K evaluated at, K=2..kmax
#' \item Y The data Y provided
#' \item kmax The maximum K requested
#' }
#' 
#' @seealso \code{\link{Clarity_fixedK}} is used at each K, and documents the important parameters. \code{\link{c_ScanOnly}} is used to generate the initial ClarityScan object independently for each K.
#' @export
#' @examples
#' \donttest{
#' scanraw=Clarity_Scan(dataraw) # Generate an initial run
#' scanraw=Clarity_Scan(dataraw,clist=scanraw) # Run it for longer
#' }
Clarity_Scan<-function(Y,
                      kmax=20,
                      clist=NULL,
                      niter=10,thresh=0.1,tmax=10000,
                      verbose=TRUE,
                      matrixdistmin=1e-8,
                      objectivedistmin=0.001,...){
    if(is.null(clist)){
        if(verbose){
            print(paste0("Creating a new clist using c_ScanOnly and kmax=",kmax,". This may take some time for larger problems."))
        }
        clist=c_ScanOnly(Y,kmax=kmax,tmax=tmax,
                               matrixdistmin=matrixdistmin,
                               objectivedistmin=objectivedistmin,
                               verbose=verbose,
                               ...)        
    }else{
        if(class(clist)!="ClarityScan") stop("clist must be of class ClarityScan as returned by Clarity_Scan")
        kmax=clist$kmax
        if(verbose){
            print(paste0("Using provided clist which uses kmax=",kmax,". "))
        }
    }
    
    if(niter>0){
        if(verbose){
            print(paste0("Performing ",niter," Backward-Forward optimisation iterations. These often start slow and speed up as the solution stabilises."))
        }
        clist=c_ForwardBackward(clist=clist,niter=niter,
                              thresh=thresh,tmax=tmax,
                              matrixdistmin=matrixdistmin,
                              objectivedistmin=objectivedistmin,
                              verbose=verbose,
                              ...)
    }
    if(verbose){
        print(paste0("If the algorithm terminated due to running out of iterations, consider calling Clarity_Scan again passing in the object it just gave you as clist."))
    }
    return(clist)
}


###############################
#' @title Predict one Similarity from a Clarity or ClarityScan object
#'
#' @description
#' Takes each Clarity object and learns a new X for the new data Ynew provided.
#' 
#' @param Ynew The data to be predicted
#' @param clist The learned Clarity object, either of class "Clarity" or of class "ClarityScan"
#' 
#' @keywords mixture
#' @return An object of the same class as clist provided, with updated Y, X and derived features.
#' 
#' @seealso \code{\link{Clarity_Scan}} to generate an estimate for clist to be used here.
#' @export
#' @examples
#' \donttest{
#' scanraw=Clarity_Scan(dataraw) # Generate an initial run
#' scanraw=Clarity_Scan(dataraw,clist=scanraw) # Run it for longer
#' # Apply it to a new dataset with the same structure
#' scanrepfromraw=Clarity_Predict(datarep,clist=scanraw) 
#' # Apply it to a new dataset with slightly different structure
#' scanmixfromraw=Clarity_Predict(datamix,clist=scanraw) 
#' }
Clarity_Predict<-function(Ynew,clist) {
    if(class(clist)=="ClarityScan") {
        ret=list()
        ret$scan=lapply(clist$scan,Clarity_Predict,Ynew=Ynew)
        ret$objectives=sapply(clist$scan,function(x)x$objective)
        ret$klist=clist$klist
        ret$Y=Ynew
        ret$kmax=clist$kmax
        class(ret)="ClarityScan"
        return(ret)
    }else if(class(clist)=="Clarity"){
        A=clist$A
        X=c_updateX(A,clist$X,Ynew)
        calXlist=c_calXfromX(A,X)
        prediction=A %*% X %*% t(A)
        objective=c_objectivefunction(A,X,Ynew)
        ret=list(Y=Ynew,A=A,X=X,
                 calX= calXlist$calX,
                 prediction = prediction, 
                 C = calXlist$C,
                 tCinv = calXlist$tCinv,
                 Afirst=A,
                 Xfirst=X,
                 calXfirst=calXlist,
                 numiters=1,
                 matrixdist=0,
                 objectivedist=objective,
                 objective=objective,
                 matrixdistmin=clist$matrixdistmin,
                 objectivedistmin=clist$objectivedistmin
                 )
        class(ret)="Clarity"
        return(ret)
    }else{
        stop("Must provide a Clarity or ClarityScan object")
    }
}


###############################
#' @title Plot a Clarity object and another predicted from it
#'
#' @description
#' Takes two Clarity objects over the same samples and plots both of predictions, their associated data, and their residuals.
#' 
#' @param c1 The first, "learned" Clarity object
#' @param c2 The second, "predicted" Clarity object
#' @param order (default NULL) the plotting order. If NULL, the data c1$Y are fit with a dendrogram to obtain a good ordering. If NA, no reordering is done. Otherwise a specified ordering can be given.
#' @param rownames (default NULL) the names to be given to the rows, in the original order. By default it is taken from rownames(c1$Y)
#' @param zlim1 (default NULL) range of the prediction and data plots for c1. Default: range of those data.
#' @param zlim2 (default NULL) as above for c1.
#' @param zlimresidual (default NULL) range of the residual plot. Default: range of residuals for c1 and c2. It is not allowed to use different scales for each residual.
#' @param name1 (default "Learned") Name for the c1 data
#' @param name2 (default "Predicted") Name for the c2 data
#' @param plotnames (default c("Fit","Data","Residuals")) names for the Prediction, Data, and Residuals plots. Set to "" to disable additional naming
#' @param cex.axis (default 1) cex for the row/column names. Set to 0 to disable axis plotting
#' @param cex.main (default 1) cex for the titles
#' @param range.cex (default 0.5) cex for the residual range (set to 0 to disable plotting)
#' @param range.line (default 0) line for the residual range text (for mtext)
#' @param ... Extra arguments for image: consider using different palettes.
#' @keywords mixture
#' @return A list containing objects from the above list of possible inputs.
#' 
#' @seealso \code{\link{Clarity_Scan}} to generate an estimate for c1 to be used here, and \code{\link{Clarity_Predict}} to generate an estimate for c2.
#' @export
#' @examples
#' \donttest{
#' scanraw=Clarity_Scan(dataraw) # Generate an initial run
#' scanraw=Clarity_Scan(dataraw,clist=scanraw) # Run it for longer
#' # Apply it to a new dataset with the same structure
#' scanrepfromraw=Clarity_Predict(datarep,clist=scanraw) 
#' # Apply it to a new dataset with slightly different structure
#' scanmixfromraw=Clarity_Predict(datamix,clist=scanraw)
#' # Plot the case where the residuals shouldn't matter
#' Clarity_ComparisonPlot(scanraw$scan[[19]],scanrepfromraw$scan[[19]],
#'       name2="Predicted (no structural change)",zlimresidual=c(-1,1))
#' # Plot the case where the residuals should matter
#' Clarity_ComparisonPlot(scanraw$scan[[19]],scanmixfromraw$scan[[19]],
#'       name2="Predicted (one structural change)",zlimresidual=c(-1,1))
#' }
Clarity_ComparisonPlot=function(c1,c2,order=NULL,rownames=NULL,
                               zlim1=NULL,zlim2=NULL,
                               zlimresidual=NULL,
                               name1="Learned",name2="Predicted",
                               plotnames=c("Fit","Data","Residuals"),
                               cex.axis=1,cex.main=1,
                               range.cex=0.5,range.line=0,...){
    ## Work out what ordering to use
    n=dim(c1$prediction)[1]
    if(any(is.null(order))){
        d1=stats::as.dendrogram(stats::hclust(stats::dist(c1$Y)))
        Rowv <- rowMeans(c1$Y)
        d1=stats::reorder(d1, Rowv)
        order=labels(d1)
    }else if(any(is.na(order))){
        order=1:n
    }
    if(any(is.null(rownames))) rownames=rownames(c1$Y)
    if(length(plotnames)<3) plotnames=rep("",3)
    
    ## Update everything with the desired order
    rownames=rownames[order]
    c1pred=(c1$prediction)[order,order]
    c2pred=(c2$prediction)[order,order]
    c1Y=(c1$Y)[order,order]
    c2Y=(c2$Y)[order,order]
    c1resid=c1pred-c1Y
    c2resid=c2pred-c2Y

    ## Sort out the scales of the plot
    if(is.null(zlim1)){
        zlim1=range(c(range(c1pred),range(c1Y)))
    }
    if(is.null(zlim2)){
        zlim2=range(c(range(c2pred),range(c2Y)))
    }
    if(is.null(zlimresidual)){
        zlimresidual=range(c(c1resid,c2resid))
    }
    realzlimresidual1=range(c1resid)
    realzlimresidual2=range(c2resid)

    ## Make the actual plots
    graphics::par(mfrow=c(2,3))
    ## Prediction for c1
    graphics::image(1:n,1:n,c1pred,xlab="",ylab="",
          axes=F,zlim=zlim1,main=paste(name1,plotnames[1]),cex.main=cex.main)
    if(cex.axis>0) graphics::axis(1,at=1:n,labels=order,las=2,cex.axis=cex.axis);
    if(cex.axis>0) graphics::axis(2,at=1:n,labels=order,las=2,cex.axis=cex.axis)

    graphics::image(1:n,1:n,c1Y,xlab="",ylab="",
         axes=F,zlim=zlim1,main=paste(name1,plotnames[2]),cex.main=cex.main)
    if(cex.axis>0)     graphics::axis(1,at=1:n,labels=order,las=2,cex.axis=cex.axis)
    if(cex.axis>0)     graphics::axis(2,at=1:n,labels=order,las=2,cex.axis=cex.axis)
    
    graphics::image(1:n,1:n,c1resid,xlab="",ylab="",
                    axes=F,zlim=zlimresidual,
                    main=paste(name1,plotnames[3]),
                    cex.main=cex.main)
    if(range.cex>0) graphics::mtext(paste0("Range (",paste(format(realzlimresidual1,digits=2),collapse=","),")"),
                                  3,cex=range.cex,line=range.line)

    if(cex.axis>0)     graphics::axis(1,at=1:n,labels=order,las=2,cex.axis=cex.axis)
    if(cex.axis>0)     graphics::axis(2,at=1:n,labels=order,las=2,cex.axis=cex.axis)

    ## Prediction for c2
    graphics::image(1:n,1:n,c2pred,xlab="",ylab="",
          axes=F,zlim=zlim2,main=paste(name2,plotnames[1]),cex.main=cex.main)
    if(cex.axis>0)     graphics::axis(1,at=1:n,labels=order,las=2,cex.axis=cex.axis)
    if(cex.axis>0)     graphics::axis(2,at=1:n,labels=order,las=2,cex.axis=cex.axis)

    graphics::image(1:n,1:n,c2Y,xlab="",ylab="",
         axes=F,zlim=zlim2,main=paste(name2,plotnames[2]),cex.main=cex.main)
    if(cex.axis>0)     graphics::axis(1,at=1:n,labels=order,las=2,cex.axis=cex.axis)
    if(cex.axis>0)     graphics::axis(2,at=1:n,labels=order,las=2,cex.axis=cex.axis)
    
    graphics::image(1:n,1:n,c2resid,xlab="",ylab="",
         axes=F,zlim=zlimresidual,main=paste(name2,plotnames[3]),cex.main=cex.main)
    if(range.cex>0) graphics::mtext(paste0("Range (",paste(format(realzlimresidual2,digits=2),collapse=","),")"),
                                  3,cex=range.cex,line=range.line)
    if(cex.axis>0)     graphics::axis(1,at=1:n,labels=order,las=2,cex.axis=cex.axis)
    if(cex.axis>0)     graphics::axis(2,at=1:n,labels=order,las=2,cex.axis=cex.axis)

    invisible(list(order=order,
                rownames=rownames,
                c1pred=c1pred,
                c2pred=c2pred,
                c1Y=c1Y,
                c2Y=c2Y,
                c1resid=c1resid,
                c2resid=c2resid,
                zlim1=zlim1,
                zlim2=zlim2,
                zlimresidual
                ))
}
