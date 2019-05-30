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



###############################
#' @title Get a plot value in a sensible way
#'
#' @description
#' Generate a default colour/line type/pch etc for a ClarityScan object. Users are suggested to supply their own choices rather than call c_getplotval directly.
#' 
#' @param x the provided value
#' @param i the index
#' @param clist the list of clarityscan objects
#' @param default what to return if x is null
#'
#' @return the desired plot symbol
c_getplotval <- function (x,i,clist,default){
    if(is.null(x))ret=default
    else if (length(x)!=length(clist)) ret=x[1]
    else ret=x[i]
    ret
}

###############################
#' @title Plot the residuals of a set of ClarityScan objects as a function of k
#'
#' @description
#' Makes a plot of (objective function) against (k) for all ClarityScan objects, automatically doing sensible things with color, line type, legend, etc. For most arguments, try passing vectors of length P=number of plots to make.
#'
#' @param clist the list of P clarityscan objects
#' @param add (default FALSE): if FALSE, call plot()
#' @param legendloc (default: "bottomleft") where the legend will be plotted. Set to non-chacter to disable legend plotting
#' @param legendbty (default: "n") box type (bty) in call to legend
#' @param name (Default: detect from clist, or use Scan1..P) names of clist items for legend
#' @param lty Line type. Default: 1
#' @param col Color. Default: 1..P
#' @param pch point type. Default: 1
#' @param lwd line width. Default: 1
#' @param lines (Default: TRUE) whether to plot lines
#' @param points (Default: TRUE) whether to plot points
#' @param xlim (Default: range of data) xlim for plot()
#' @param ylim (Default: range of data) ylim for plot()
#' @param xlab Label for x axis
#' @param ylab Label for y axis
#' @param main Label for main
#' @param extern (Default: TRUE) whether being called externally or not.
#' @param ... Additional parameters to plot
#' @return the desired plot symbol
#' @export
#' 
Clarity_ObjectivePlot <-function(clist,add=FALSE,
                                 legendloc="bottomleft",
                                 legendbty="n",
                                 name=NULL,
                                 lty=NULL,col=NULL,pch=NULL,lwd=NULL,
                                 lines=TRUE,points=TRUE,
                                 xlim=NULL,ylim=NULL,
                                 xlab="K",
                                 ylab="Objective Function",
                                 main="",
                                 extern=TRUE,
                                 ...){
    if(class(clist)%in%c("ClarityScan","ClarityScanExtend")) {
        if(extern){
            print("Warning: called Clarity_ObjectivePlot on a ClarityScan object. Expected a list of ClarityScan objects with names. Will proceed with the name 'scan'.");
            clist=list("scan"=clist)
        }else{
            if(!add){
                if(is.null(xlim)) xlim=range(clist$klist)
                if(is.null(ylim)) xlim=range(clist$objectives)
                if(is.null(lty)) lty=1
                if(is.null(col)) col=1
                if(is.null(pch)) pch=1
                if(is.null(lwd)) lwd=1
                graphics::plot(clist$klist,clist$objectives,xlim=xlim,ylim=ylim,
                               type="n",xlab=xlab,ylab=ylab,main=main,...)
            }
            if(lines)graphics::lines(clist$klist,clist$objectives,lty=lty,lwd=lwd,col=col)
            if(points)graphics::points(clist$klist,clist$objectives,pch=pch,col=col)
            ret=list(name=name,lty=lty,lwd=lwd,col=col,pch=pch)
            return(ret)
        }
    }
    if(class(clist)=="list"){
        if(class(clist[[1]])%in%c("ClarityScan","ClarityScanExtend")){
            if(is.null(xlim)) xlim=range(sapply(clist,function(x)x$klist))
            if(is.null(ylim)) ylim=range(sapply(clist,function(x)x$objectives))
            if(is.null(name) && !is.null(names(clist))) name=names(clist)
            Clarity_ObjectivePlot(clist[[1]],add=add,lines=FALSE,points=FALSE,
                                           xlim=xlim,ylim=ylim,extern=FALSE,...)
            if(length(clist)>=1) {
                for(i in 1:length(clist)){
                    tcol=c_getplotval(col,i,clist,i)
                    tlwd=c_getplotval(lwd,i,clist,1)
                    tlty=c_getplotval(lty,i,clist,1)
                    tpch=c_getplotval(pch,i,clist,1)
                    tname=c_getplotval(name,i,clist,paste0("Scan",i))
                    if(i==1) ret=list()
                    ret[[i]]=Clarity_ObjectivePlot(clist[[i]],add=TRUE,
                                                    name=tname,
                                                    lines=lines,points=points,
                                                    col=tcol,lwd=tlwd,lty=tlty,pch=tpch,extern=FALSE)
                }
            }
        }else stop("Unrecognised structure for clist")
        if(class(legendloc)=="character")
            graphics::legend(legendloc,legend=sapply(ret,function(x)x$name),
               lty=sapply(ret,function(x)x$lty),
               col=sapply(ret,function(x)x$col),
               text.col=sapply(ret,function(x)x$col),
               lwd=sapply(ret,function(x)x$lwd),
               pch=sapply(ret,function(x)x$pch),
               bty=legendbty)

        return(invisible(ret))
    }else{
        stop("Unrecognised structure for clist")
    }
}
