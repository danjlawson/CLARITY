
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
