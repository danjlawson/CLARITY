
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
#' @seealso \code{\link{c_Bootstrap}}, which uses this code.
c_Procrust=function(Ystart,Ytarget){
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
    ret[ret<0]=0
    colnames(ret)=colnames(Ystart)
    ret
}


###############################
#' @title Perform a Mean Scale transformation of one matrix into another
#'
#' @description
#'
#' Scale the mean of the matrix provided in the first argument, into the target provided in the second.
#' @param Ystart the matrix to transform
#' @param Ytarget The target matrix
#' @return A matrix which is the best fitting transformation of Ystart to Ytarget
#' @export
#' @seealso \code{\link{c_Bootstrap}}, which uses this code.
c_MeanScale=function(Ystart,Ytarget){
    if (nrow(Ytarget) != nrow(Ystart)) {
        stop("Dimensions differ between Ystart and Ytarget")
    }
    if (ncol(Ytarget) != ncol(Ystart)) {
        stop("Dimensions differ between Ystart and Ytarget")
    }
    Ystart*mean(Ytarget)/mean(Ystart)
}

###############################
#' @title Identity function taking two arguments
#'
#' @description
#'
#' Returns the Y matrix untransformed. Provided for input as the transform variable in \code{\link{c_Bootstrap}}.
#' @param Ystart the matrix to transform
#' @param Ytarget The target matrix
#' @return Ystart, unmodified
#' @export
#' @seealso \code{\link{c_Bootstrap}}, which uses this code.
c_Identity=function(Ystart,Ytarget=NULL){
    Ystart
}

