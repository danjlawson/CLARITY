
###############################
#' @title Merge a rectangular matrix into a lower dimensional matrix using mixtures
#'
#' @description
#' Takes a matrix x (N by P) and a mixture matrix A (N by K) and computes an N by K matrix:
#' \code{t(t(x) \%*\% A)}
#' Also has the facility to avoid the final transpose, and to replace A by normalised version that switches the computation from a weighted sum across "clusters" to a weighted mean by using
#' \code{An=t(t(A)/colSums(A))}
#' instead of A.
#'
#' @param x a matrix of dimension N by P
#' @param A a mixture matrix  of dimension N by K
#' @param rotate (default=FALSE) whether to return a rotated version of the result (by omitting the final transpose)
#' @param type (default="mean") whether to reaturn the "mean" or the "sum" of pairwise elements across the (soft) clustering described by A
#'
#' @return A matrix of dimension N by K
#' @seealso \code{\link{plot.ClarityScan}} usez this to make a summary of Clarity persistences.
#' @export
c_Merge=function(x,A,rotate=FALSE,type="mean"){
    if(type=="mean"){
        An=t(t(A)/colSums(A))
        ret=(t(x) %*% An)
    }else if(type=="sum"){ 
        ret=(t(x) %*% A)
    }else{
        stop("type must be either \"mean\" or \"sum\"")
    }
    if(rotate==FALSE){
        ## We rotated for the calculation, so rotate back
        ret=t(ret)
    }
    return(ret)
}

###############################
#' @title Merge a square matrix into a square lower dimensional square matrix using mixtures
#'
#' @description
#' Takes a matrix x [N by N] and a mixture matrix A (N by K) and computes a K by K matrix:
#' \code{t(t(A) \%*\% t(x) \%*\% A)}
#' Also has the facility to avoid the final transpose, and to replace A by normalised version that switches the computation from a weighted sum across "clusters" to a weighted mean by using
#' \code{An=t(t(A)/colSums(A))}
#' instead of A.
#'
#' @param x a matrix of dimension N by N
#' @param A a mixture matrix  of dimension N by K
#' @param rotate (default=FALSE) whether to return a rotated version of the result (by omitting the final transpose)
#' @param type (default="mean") whether to reaturn the "mean" or the "sum" of pairwise elements across the (soft) clustering described by A
#'
#' @return A matrix of dimension K by K
#' @seealso \code{\link{plot.Clarity}} usez this to make a summary of Clarity residuals.
#' @export
c_MergeSquare=function(x,A,rotate=FALSE,type="mean"){
    if(type=="mean"){
        An=t(t(A)/colSums(A))
        ret=(t(An) %*% t(x) %*% An)
    }else if(type=="sum"){ 
        ret=(t(A) %*% t(x) %*% A)
    }else{
        stop("type must be either \"mean\" or \"sum\"")
    }
    if(rotate==FALSE){
        ## We rotated for the calculation, so rotate back
        ret=t(ret)
    }
    return(ret)
}

###############################
#' @title Random dirichlet RV
#'
#' @description
#' Simulate n Dirichlet Random Variables with parameter alpha
#' 
#' @param n number of random dirichlet variables
#' @param alpha dirichlet parameter
#' 
#' @return A matrix with each random variate in a row
#' 
#' @export
#' 
rdirichlet<-function (n, alpha) 
{
    l <- length(alpha)
    x <- matrix(stats::rgamma(l * n, alpha), ncol = l, byrow = TRUE)
    sm <- x %*% rep(1, l)
    x/as.vector(sm)
}
###############################
#' @title Match arguments
#'
#' @description
#' find the best partial match for x in y, ignoring case
#' 
#' @param x the argument to match fuzzily
#' @param y a character vector containing the options
#' 
#' @return The index of the matching value (NA if no valid match)
#' 
#' @export
#' 
c_argmatch=function(x,y){
    ret=pmatch(tolower(x),tolower(y),nomatch=0)
    if(ret>0) return(y[ret])
    "NA"
}
