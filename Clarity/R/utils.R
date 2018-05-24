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
