###############################
#' @title Get a centering matrix of specified dimension 
#'
#' @description
#' Gets a centering matrix Cn so that Cn %*% Y is centered.
#' @param n either the dimension of Y, or Y itself.
#' @return A matrix with same dimensions as Y, centered.
#' @export
c_getC=function(n){ # The centering matrix from
    if(is(n,"matrix")) n=dim(n)[1]
    diag(n) - (1/n)*matrix(1,nrow=n,ncol=n)
}
###############################
#' @title (Double) Center a matrix
#'
#' @description
#' Centers a matrix by rows and/or columns, nominally computing 'Cn %*% Y %*% Cm' for an n by m matrix Y
#'
#' @param Y the matrix to be centered
#' @param col (default: TRUE) Whether to center the columns to ensure they have mean zero
#' @param row (default: TRUE) Whether to center the rows to ensure they have mean zero
#' @param method Either "mean" (Default) or "matrix"; Whether to center with the theoretically nice centering matrix \code{\link{c_getC}} or just subtracting row and column means appropriately, which is quicker, and supports na values.
#' @param na.rm Whether to remove NA values if method=="mean"
#' @return A matrix with same dimensions as Y, centered.
#' @export
c_Center <- function (Y, col = TRUE, row = TRUE, method = "mean",na.rm=TRUE) 
{
    ret = Y
    if (col) {
        if (method == "matrix") {
            Cn1 = c_getC(dim(Y)[1])
            ret = Cn1 %*% Y
        }
        else if (method == "mean") {
            ret = ret - rowMeans(ret,na.rm=na.rm)
        }
        else stop("Invalid method: must be \"matrix\" or \"mean\"")
    }
    if (row) {
        if (method == "matrix") {
            Cn2 = c_getC(dim(Y)[2])
            ret = ret %*% Cn2
        }
        else if (method == "mean") {
            ret = t(t(ret) - colMeans(ret,na.rm=na.rm))
        }
        else stop("Invalid method: must be \"matrix\" or \"mean\"")
    }
    rownames(ret) = rownames(Y)
    colnames(ret) = colnames(Y)
    ret
}
###############################
#' @title Perform a Procrustes transformation of one matrix into another
#'
#' @description
#'
#' Procrust the matrix provided in the first argument, into the target provided in the second.  Performs Translation and Dilation.
#' @param Ystart the matrix to transform
#' @param Ytarget The target matrix
#' @param positive (default=TRUE) whether to force any negative components to zero. This is the correct default for similarities, but not for other uses.
#' @return A matrix which is the best fitting transformation of Ystart to Ytarget
#' @export
#' @seealso \code{\link{c_Bootstrap}}, which uses this code.
c_Procrust<- function (Ystart, Ytarget,positive=TRUE) 
{
    if (nrow(Ytarget) != nrow(Ystart)) {
        stop("Dimensions differ between Ystart and Ytarget")
    }
    if (ncol(Ytarget) != ncol(Ystart)) {
        stop("Dimensions differ between Ystart and Ytarget")
    }
    n <- nrow(Ystart)
    m <- ncol(Ystart)
    J <- diag(n) - 1/n * matrix(1, n, n)
    C <- t(Ytarget) %*% J %*% Ystart
    svd.out <- svd(C)
    R <- svd.out$v %*% t(svd.out$u)
    s <- 1
    mat1 <- t(Ytarget) %*% J %*% Ystart %*% R
    mat2 <- t(Ystart) %*% J %*% Ystart
    s.numer <- 0
    s.denom <- 0
    for (i in 1:m) {
        s.numer <- s.numer + mat1[i, i]
        s.denom <- s.denom + mat2[i, i]
    }
    s <- s.numer/s.denom
    tt <- matrix(0, m, 1)
    tt <- 1/n * t(Ytarget - s * Ystart %*% R) %*% matrix(1, n, 1)
    ret = s * Ystart %*% R + matrix(tt, n, m, byrow = TRUE)
    if(positive) ret[ret < 0] = 0
    colnames(ret) = colnames(Ystart)
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
    tmp=Ystart-mean(Ystart)
    tmp=tmp*stats::sd(Ytarget)/stats::sd(Ystart)
    tmp=tmp+mean(Ytarget)
    tmp
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

