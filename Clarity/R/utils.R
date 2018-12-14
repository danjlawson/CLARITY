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
