###############################
#' @title Plot a matrix in the Clarity "Persistence Chart" style
#'
#' @description
#' Takes a matrix and plots it as an advanced heatmap, including markup for significance, explicit scaling of colours, text on elements, and more.
#'
#' NB The recommended interface is  \code{\link{plot.Clarity}} and \code{\link{plot.ClarityScan}}.
#'
#' @param x The matrix object to be drawn
#' @param signif (default=NULL) an optional logical matrix of the same dimension as x. If provided, instead of a regular image, entries that are FALSE are drawn smaller and faded. If not provided, the plot is drawn as if entries are TRUE.
#' @param mar (default=c(6,6,4,1) margins of the plot as in \code{\link{par}}(mar). Set to NULL to disable setting them.
#' @param scalefun (default=sqrt) scaling function for the entries. The default is appropriate when plotting squared values such as  square residuals, persistences, etc. function(x){log(1+x)} can often be a good choice, as can the identity I.
#' @param zlim (default=NULL) range of the colour scale. Defaults to the range of x. If zlim=="symmetric" then it is symmetric around 0.
#' @param cols (default=colorRampPalette(c("white","yellow","orange","red","darkred"))(100)) a colour scale. It is essential if signif is provided that the colors are provided in "#RRGGBBAA" for red/green/blue/alpha hex values.
#' @param text (default=FALSE) whether to write the values on each matrix element. This is useful for small matrices.
#' @param cex.text (default=1) text size if text==TRUE.
#' @param digits (default=2) number of digits to write matrix values to if text==TRUE.
#' @param axes (default=TRUE) whether to draw any axes.
#' @param las.axis (default=1) \code{\link{par}}(las) for axis text.
#' @param cex.axis (default=1) \code{\link{par}}(cex) for axis text.
#' @param cex.axis.X (default=NULL) cex for X axis. Inherits from cex.axis if NULL. Set to 0 to disable X axis.
#' @param cex.axis.Y (default=NULL) cex for Y axis. Inherits from cex.axis if NULL. Set to 0 to disable Y axis.
#' @param axis.top (default=FALSE) whether to draw the X axis at the top, rather than the bottom, of the image.
#' @param etext (default=NULL) an optional matrix of additional text to be drawn in each cell underneath the main text.
#' @param etextgap (default=0.25) distance below the center of each rectangle to place the etext content.
#' @param cex.etext (default=1) size of additional text.
#' @param userect (default=FALSE) whether to use the "rectangle" plotting style, or the image plotting style. This is overriddent to TRUE if signif is provided.
#' @param rectdelta (default=c(0.8,0.8)) scaling factor for (X,Y) for non-significant cells.
#' @param signiffade (default="55") fading in the alpha channel for non-significant cells.
#' @param imageback (default="white") for the rectangle plotting style, the background colour that will appear around the outside of non-significant cells.
#' @param tol (default=1e-10) for the rectangle plotting style, the tolerance for determining which colour bin entries should appear
#' @param main (default="") Title as in \code{\link{image}}.
#' @param xlab (default="") Title as in \code{\link{image}}.
#' @param ylab (default="") Title as in \code{\link{image}}.
#' @param ... Additional parameters to \code{\link{image}}
#'
#' @return NULL, invisibly. This funciton is used for its side effects in plotting.
#'
#' @seealso \code{\link{plot.Clarity}}, \code{\link{plot.ClarityScan}} require much less manual effort to extract various matrices from Clarity and ClarityScan objects. \code{\link{c_legend}} for one way to present a legend.
#' @examples
#' \donttest{
#' scan=Clarity_Scan(dataraw)
#' predmix=Clarity_Predict(datamix,scan)
#' scanpred=Clarity_Predict(datamix,scan)
#' scanbootstrap=Clarity_Bootstrap(scan,target=datamix,D=datarawD)
#'
#' ## The recommended way to plot:
#' plot(scanpred,signif=scanbootstrap)
#'
#' ## The manual way to plot:
#' ## Extract observed persistences
#' P=Clarity_Persistence(predmix)
#' ## Compute pvalues
#' pvals=Clarity_BScompare(scanbootstrap,P)
#' ## pvals is a matrix of dimension N by K
#' signif=pvals<0.01
#'
#' ## Now pass to Clarity_Chart
#' Clarity_Chart(P,signif)
#' }
#' @export
#' 
Clarity_Chart<-function(x,
                        signif=NULL,
                        mar=c(6,6,4,1),
                        scalefun=sqrt,
                        zlim=NULL,
                        cols=grDevices::colorRampPalette(c("white","yellow","orange","red","darkred"))(100),
                        text=FALSE,
                        cex.text=1,
                        digits=2,
                        axes=TRUE,
                        las.axis=1,
                        cex.axis=1,
                        cex.axis.X=NULL,cex.axis.Y=NULL,
                        axis.top=FALSE,
                        etext=NULL,
                        etextgap=0.25,
                        cex.etext=1,
                        userect=FALSE,
                        rectdelta=c(0.8,0.8),
                        signiffade="55",
                        imageback="white",
                        tol=1e-10,
                        main="",
                        xlab="",ylab="",...){
    if(!any(is.null(mar))) graphics::par(mar=mar)
    if(length(rectdelta)==1)rectdelta=rep(rectdelta,2)
    if(is.null(zlim)) zlim=range(stats::na.omit(as.numeric(scalefun(x))))
    if(class(zlim)=="character" && zlim=="symmetric"){
        zlim=max(abs(stats::na.omit(as.numeric(scalefun(x)))))
        zlim=c(-zlim,zlim)
    }
    if(any(!is.null(signif))) {
        if(class(signif)=="numeric"){
            signif=matrix(TRUE,nrow=dim(x)[1],ncol=dim(x)[2])
        }
        userect=TRUE
    }
    if(userect){
        imagecols=imageback
        if(all(is.null(signif))){
            signif=x
            signif[]=TRUE
        }
    }else{
        imagecols=cols
    }
    graphics::image(1:(dim(x)[1]+1)-0.5,
          1:(dim(x)[2]+1)-0.5,
          scalefun(x),col=imagecols,
          axes=F,xlab=xlab,ylab=ylab,
          zlim=zlim,main=main,...)
    if(userect) {
        tthresh=seq(zlim[1]-tol,zlim[2]+tol,length.out=1+length(cols))
        rectcols=x
        rectcols[]=sapply(x,function(x1){
            if(is.na(x1)) return(0)
            min(which(tthresh>=scalefun(x1)))
        })-1
        rcols=paste0(substr(cols,1,7),signiffade)
        for(i in 1:dim(x)[1])
            for(j in 1:dim(x)[2]){
                tcol=ifelse(signif[i,j],cols[rectcols[i,j]],rcols[rectcols[i,j]])
                graphics::rect(i-0.5*ifelse(signif[i,j],1,rectdelta[1]),
                    j-0.5*ifelse(signif[i,j],1,rectdelta[2]),
                    i+0.5*ifelse(signif[i,j],1,rectdelta[1]),
                    j+0.5*ifelse(signif[i,j],1,rectdelta[2]),
                    col=tcol,border=NA)
            }
    }
    if(axes){
        if(is.null(cex.axis.X))cex.axis.X=cex.axis
        if(is.null(cex.axis.Y))cex.axis.Y=cex.axis
        if(cex.axis.X>0){
            if(axis.top){taxis=3}else{taxis=1}
            graphics::axis(taxis,at=1:dim(x)[1],
                           labels=rownames(x),las=las.axis,
                           cex.axis=cex.axis.X)
        }
        if(cex.axis.Y>0){
            graphics::axis(2,at=1:dim(x)[2],
                           labels=colnames(x),las=las.axis,
                           cex.axis=cex.axis.Y)
        }
    }
    if(text){
        for(i in 1:dim(x)[1]) {
            graphics::text(i,1:dim(x)[2],
                 round(x,digits=digits)[i,],
                 cex=cex.text)
        }
    }
    if(!is.null(etext)) if(class(etext)=="matrix") if(all(dim(etext)==dim(x)))
    {
        for(i in 1:dim(x)[1]) {
            text(i,1:dim(x)[2]-etextgap,etext[i,],cex=cex.etext)
        }
    }
    invisible(NULL)
}


###############################
#' @title Plot residuals from a Clarity object
#'
#' @description
#' Plot a Residuals Chart, which is an image representing which subjects are poorly explained by a Clarity fit.
#' 
#' This is the main plot fuction for residuals. It has several convenience manipulations to extract residuals and statistical significance associated with them via bootstrapping, as well as summarising them by clustering or mixtures.
#'
#' @param x A Clarity object.
#' @param signif (default=NULL) either a ClarityBootstrap object as returned by \code{\link{Clarity_Bootstrap}}; or a matrix of p-values as returned by \code{\link{Clarity_BScompare}}; or a logical matrix of elements that are significant (the same dimensions as the similarities modelled in x). In all cases it is transformed to a logical matrix and passed to \code{\link{Clarity_Chart}}.
#' @param A (default=NULL) A mixture matrix by which to summarise the matrix using \code{\link{c_Merge}}. This might e.g. represent a clustering of the elements in the similarities. It must be an N (number of subjects) by K' (number of desired clusters) matrix with named columns.
#' @param order (default=NULL) preferred order. If NULL, use the provided ordering. Otherwise, either a numeric vector of length N, or a permutation of the rownames of the original similarity Y.
#' @param plot (default=TRUE) whether to actually create the plot or store the computed matrices for later plotting with \code{\link{plot.ClarityPlot}}.
#' @param rotate (default=FALSE) whether to rotate the Residuals.
#' @param type (default="mean") if A is not null, should we report the "mean" or the "sum" of the elements in each pair of (soft) clusters?
#' @param thresh (default=0.01) threshold for translating bootstrapped p-values into significant/not-significant calls for the matrix signif.
#' @param diag (default=NA) how to treat the diagonal of residuals when plotting residuals, as used by \code{\link{Clarity_Extract}}. Set to NULL to use the observed values, NA to treat them as missing, or 0.
#' @param summary (default=abs) how to summarise residuals for plotting (not p-values), as used by \code{\link{Clarity_Extract}}. Set to I to extract the actual residuals.
#' @param population (default="bestcol") how to define the null model when computing p-values from bootstrapped data using \code{\link{Clarity_BScompare}}. The default is conservative by comparing each element to the most extreme values in its row, but "element" can be used to compare each element to only bootstrapped versions of itself.
#' @param verbose (default=FALSE) whether to output information about what transformations are being done.
#' @param ... Additional parameters to \code{\link{Clarity_Chart}}.
#'
#' @return A ClarityPlot object (invisibly) which is a list containing:
#' \itemize{
#' \item clist The provided ClarityScan object
#' \item A The provided A matrix
#' \item R The extracted Persistence matrix, transformed and merged
#' \item P0 The raw Persistence matrix as returned from \code{\link{Clarity_Persistence}}
#' \item signif0 The raw logical significance matrix, either as provided or as returned from \code{\link{Clarity_BScompare}} if provided with a ClarityScanBootstrap object (same dimension as P0)
#' \item signif The significance matrix, transformed and merged (same dimension as P)
#' }
#' @seealso \code{\link{Clarity_Chart}} which ultimately does the plotting; \code{\link{plot.ClarityPlot}} which is called to perform that plotting; \code{\link{c_MergeSquare}} which is used to perform simplification using the mixture matrix A. \code{\link{c_legend}} for one way to present a legend.
#' @examples
#' \donttest{
#' scan=Clarity_Scan(dataraw) ## Core Clarity
#' predmix=Clarity_Predict(datamix,scan) ## Core prediction
#' ## Bootstrap residuals
#' k10bootstrap=Clarity_Bootstrap(Clarity_Extract(scan,10),
#'                                target=datamix, 
#'                                D=datarawD)
#' ## Plot Residuals
#' plot(Clarity_Extract(predmix,10),
#'               signif=k10bootstrap)
#' ## Plot merged clusters
#' plot(Clarity_Extract(predmix,10),
#'               signif=k10bootstrap,A=datarawA)
#' 
#' }
#' @export
plot.Clarity=function(x,
                      signif=NULL,
                      A=NULL,
                      order=NULL,
                      plot=TRUE,
                      rotate=FALSE,
                      thresh=0.01,
                      type="mean",
                      diag=NA,
                      summary=abs,
                      population="bestcol",
                      verbose=FALSE,
                      ...){
    if(!any(is.null(signif))) {
        if(class(signif)=="ClarityBootstrap"){
            if(verbose)print("Extracting Residuals from provided ClarityBootstrap")
            R0=Clarity_Extract(x,diag=0)
            signif=Clarity_BScompare(signif,R0,population=population)
        }
        if(class(signif)=="matrix"){
         if(class(signif[1,1])=="numeric") {
                if(verbose)print(paste("Applying significance threshold",thresh))
                signif=signif<thresh
         }else if(verbose)print("Provided with significance matrix")
        }else stop("Invalid class of signif object, must either be a ClarityBootstrap or a matrix")
    }
    signif0=signif
    if(!any(is.null(A))) {
        if(verbose)print("Extracting residuals")
        R0=Clarity_Extract(x,diag=0)
        if(verbose)print("Merging residuals with provided with mixture matrix A")
        R=c_MergeSquare(R0,
                        A,rotate=rotate,type=type)
        if(!any(is.null(signif))) {
            if(verbose)print("Merging significance with provided with mixture matrix A")
            signif=c_MergeSquare(signif,A,rotate=rotate,type=type)
        }
    }else{
        if(verbose)print("Extracting residuals")
        R0=R=Clarity_Extract(x,diag=diag,summary=summary)
    }
    ret=list(clist=x,
             A=A,
             R=R,
             R0=R0,
             signif0=signif0,
             signif=signif,
             rotate=rotate,
             order=order)
    class(ret)="ClarityPlot"
    if(plot) {
        if(verbose)print("Plotting...")
        plot(ret,...)
    }
    return(invisible(ret))
}

###############################
#' @title Plot persistences from a ClarityScan object
#'
#' @description
#' Plot a Persistence Chart, which is an image representing which subjects are poorly explained by a Clarity fit.
#' 
#' This is the main plot fuction for persistences. It has several convenience manipulations to extract persistences and statistical significance associated with them via bootstrapping, as well as summarising them by clustering or mixtures.
#'
#' @param x A ClarityScan object.
#' @param signif (default=NULL) either a ClarityScanBootstrap object as returned by \code{\link{Clarity_Bootstrap}}; or a matrix of p-values as returned by \code{\link{Clarity_BScompare}}; or a logical matrix of elements that are significant (the same dimensions as the similarities modelled in x). In all cases it is transformed to a logical matrix and passed to \code{\link{Clarity_Chart}}.
#' @param A (default=NULL) A mixture matrix by which to summarise the matrix using \code{\link{c_MergeSquare}}. This might e.g. represent a clustering of the elements in the similarities. It must be an N (number of subjects) by K' (number of desired clusters) matrix with named columns.
#' @param order (default=NULL) preferred order. If NULL, use the provided ordering. Otherwise, either a numeric vector of length N, or a permutation of the rownames of the original similarity Y.
#' @param rotate (default=FALSE) whether to rotate the Persistences. If FALSE (default) then the complexity K is the Y-axis. If TRUE then the complexity K is the X-axis.
#' @param plot (default=TRUE) whether to actually create the plot or store the computed matrices for later plotting with \code{\link{plot.ClarityPlot}}.
#' @param type (default="mean") if A is not null, should we report the "mean" or the "sum" of the elements in each pair of (soft) clusters?
#' @param kmax (default=NULL) if provided, persistences above kmax are not included in the plot.
#' @param thresh (default=0.01) threshold for translating bootstrapped p-values into significant/not-significant calls for the matrix signif.
#' @param diag (default=0) how to treat the diagonal of residuals when computing persistences, as used by \code{\link{Clarity_Persistence}}. Set to NULL to use the observed values, but the default 0 is recommended.
#' @param f (default="RowSumsSquared") how to summarise residuals into persistences, as used by \code{\link{Clarity_Persistence}}.
#' @param what (default="Yresid") what to extract from each complexity when computing persistences, as used by \code{\link{Clarity_Persistence}}.
#' @param summary (default=abs) how to summarise residuals into persistences, as used by \code{\link{Clarity_Persistence}}.
#' @param population (default="bestcol") how to define the null model when computing p-values from bootstrapped data using \code{\link{Clarity_BScompare}}. The default is conservative by comparing each element to the most extreme values in its row, but "element" can be used to compare each element to only bootstrapped versions of itself.
#' @param verbose (default=FALSE) whether to output information about what transformations are being done.
#' @param ... Additional parameters to  \code{\link{Clarity_Chart}}.
#'
#' @return A ClarityPlot object (invisibly) which is a list containing:
#' \itemize{
#' \item clist The provided ClarityScan object
#' \item A The provided A matrix
#' \item P The extracted Persistence matrix, transformed and merged
#' \item P0 The raw Persistence matrix as returned from \code{\link{Clarity_Persistence}}
#' \item signif0 The raw logical significance matrix, either as provided or as returned from \code{\link{Clarity_BScompare}} if provided with a ClarityScanBootstrap object (same dimension as P0)
#' \item signif The significance matrix, transformed and merged (same dimension as P)
#' \item rotate logical indicating whether we rotated the matrices
#' }
#' @seealso Uses \code{\link{Clarity_Chart}} for plotting. \code{\link{plot.ClarityPlot}} is called to perform that plotting; \code{\link{c_Merge}} which is used to perform simplification using the mixture matrix A. \code{\link{c_legend}} for one way to present a legend. \code{\link{Clarity_ObjectivePlot}} for plotting the objective function.
#' @examples
#' \donttest{
#' scan=Clarity_Scan(dataraw) ## Core Clarity
#' predmix=Clarity_Predict(datamix,scan) ## Core prediction
#' scanbootstrap=Clarity_Bootstrap(scan,
#'                                 target=datamix,
#'                                 D=datarawD)
#' 
#' ## Plotting:
#' plot(predmix,signif=scanbootstrap)
#' ## Merged clusters:
#' plot(predmix,signif=scanbootstrap,A=datarawA)
#' }
#' @export
plot.ClarityScan=function(x,
                          signif=NULL,
                          A=NULL,
                          order=NULL,
                          rotate=FALSE,
                          plot=TRUE,
                          type="mean",
                          kmax=NULL,
                          thresh=0.01,
                          diag=0,
                          f = "RowSumsSquared",
                          what = "Yresid",
                          summary = abs,
                          population="bestcol",
                          verbose=FALSE,
                          ...){
    if(!any(is.null(signif))) {
        if(class(signif)=="ClarityScanBootstrap"){
            if(verbose)print("Extracting Persistence from provided ClarityScanBootstrap")
            P0=Clarity_Persistence(x,f=f,what=what,summary=abs,diag=diag)
            signif=Clarity_BScompare(signif,P0,population=population)
        }
        if(class(signif)=="matrix"){
         if(class(signif[1,1])=="numeric") {
                if(verbose)print(paste("Applying significance threshold",thresh))
                signif=signif<thresh
         }else if(verbose)print("Provided with significance matrix")
        }else stop("Invalid class of signif object, must either be a ClarityScanBootstrap or a matrix")
    }
    signif0=signif
    if(!any(is.null(A))) {
        if(verbose)print("Computing Persistence")
        P0=Clarity_Persistence(x,f=f,what=what,summary=abs,diag=0)
        if(verbose)print("Merging persistence with provided with mixture matrix A")
        P=c_Merge(P0,A,rotate=rotate,type=type)
        if(!any(is.null(signif))) {
            if(verbose)print("Merging significance with provided with mixture matrix A")
            signif=c_Merge(signif,A,rotate=rotate,type=type)
        }
    }else{
        if(verbose)print("Computing Persistence")
        P0=P=Clarity_Persistence(x,f=f,what=what,summary=abs,diag=diag)
        if(rotate) {
            if(verbose)print("Rotating persistence")
            P0=P=t(P)
            if(!any(is.null(signif))) {
                if(verbose)print("Rotating significance")
                signif=t(signif)
            }
        }
    }
    if(!is.null(kmax)){
        if(verbose)print(paste("Restricting to top",kmax,"dimensions"))
        if(rotate){
            P=P[1:kmax,,drop=FALSE]
            if(!any(is.null(signif))){
                signif=signif[1:kmax,,drop=FALSE]
            }
        }else{
            P=P[,1:kmax,drop=FALSE]
            if(!any(is.null(signif))){
                signif=signif[,1:kmax,drop=FALSE]
            }
        }
    }
    ret=list(clist=x,
             A=A,
             P=P,
             P0=P0,
             signif0=signif0,
             signif=signif,
             rotate=rotate,
             order=order)
    class(ret)="ClarityPlot"
    if(plot) {
        if(verbose)print("Plotting...")
        plot(ret,...)
    }
    return(invisible(ret))
}

###############################
#' @title Plot a ClarityPlot object
#'
#' @description
#' This is a convenience wrapper to \code{\link{Clarity_Chart}} for ClarityPlot objects returned from \code{\link{plot.Clarity}} or \code{\link{plot.ClarityScan}}. It will extract the Persistence P or Squared Residuals R, along with any signif matrix.
#'
#' @param x A ClarityPlot object
#' @param order (default=NULL) preferred order. If NULL, use the provided ordering. Otherwise, either a numeric vector of length N, or a permutation of the rownames of the original similarity Y.
#' @param ... Extra parameters to be passed to \code{\link{Clarity_Chart}}
#' 
#' @return NULL, invisibly. This funciton is used for its side effects in plotting.
#' @seealso \code{\link{Clarity_Chart}} for the actual plotting, \code{\link{plot.Clarity}} for plotting Residuals, and \code{\link{plot.ClarityScan}} for plotting Persistences. \code{\link{c_legend}} for one way to present a legend.
#' 
#' @examples
#' \donttest{
#' scan=Clarity_Scan(dataraw) ## Core Clarity
#' predmix=Clarity_Predict(datamix,scan) ## Core prediction
#' scanbootstrap=Clarity_Bootstrap(scan,target=datamix,D=datarawD)
#' 
#' ## Plotting: generate a plotting object
#' predplot=plot(predmix,signif=scanbootstrap,plot=FALSE)
#' ## Plot it
#' plot(predplot)
#' }
#' @export
#' 
plot.ClarityPlot=function(x,order=NULL,...){
    if(is.null(order))order=x$order
    if(any(!is.null(x$P))){
        val=x$P
        if(any(!is.null(order))){
            if(x$rotate){
                val=val[,order]
                x$signif=x$signif[,order]
            }else{
                val=val[order,]
                x$signif=x$signif[order,]
            }
        }
    }else if(any(!is.null(x$R))){
        val=x$R
        if(any(!is.null(order))){
            val=val[order,order]
        }
    }else{
        stop("Cannot find P or R in ClarityPlot")
    }
    invisible(Clarity_Chart(val,signif=x$signif,...))
}


###############################
#' @title Plot a legend that describes how CLARITY shows significance
#'
#' @description
#' This is a rather ugly, manual hack to add a legend outside of the main image area. It puts example rectangles that are significant and non-significant at a given location in the image co-ordinates.
#' 
#' You might want to show this some other way, but this can work with manual tuning of x,y, and size in image units.
#' 
#' @param x the x location of the legend (center of the top element); you can specifiy top and bottom elements locations by giving two values.
#' @param y the y location of the legend (center of the top element); you can specifiy top and bottom elements locations by giving two values.
#' @param size (default=c(0.5,0.5)) the size of the element as a distance from the center specified as (x, y)
#' @param rectdelta (default=c(0.8,0.8)) the relative reduction of the non-significant rectangle
#' @param imageback (default="white") colour to show behind the non-significant rectangle
#' @param border (default="grey") colour of the border to place around each rectangle. Set to NA to omit.
#' @param cex.text (default=1.5) size of the legend text
#' @param gap (default=0) gap between the two elements, if their locations were not both specified
#' @param laboffset (default=0.1) distance between the edge of the rectangle and the start of the text
#' @param col (default="#FF0000FF") base colour for significance
#' @param signiffade (default="55") change to the alpha channel for non-significance
#' @param xpd (default=NA) how to protect plotting outside of the region; see \code{\link{par}}(xpd).
#' @param text (default=c("Significant","Not Significant")) what to write as the legends.
#' @examples
#' \donttest{
#' scan=Clarity_Scan(dataraw) ## Core Clarity
#' predmix=Clarity_Predict(datamix,scan) ## Core prediction
#'
#' ## Plot bootstrapped persistence chart
#' scanbootstrap=Clarity_Bootstrap(scan,target=datamix,D=datarawD)
#' predplot=plot(predmix,signif=scanbootstrap)
#' ## Add the legend outside of the regular image space
#' c_legend(0,22,size=c(10,0.5))
#' }
#' @export
c_legend=function(x,
                  y,
                  size=c(0.5,0.5),
                  rectdelta=c(0.8,0.8),
                  imageback="white",border="grey",
                  cex.text=1.5,
                  gap=0,
                  laboffset=0.1,
                  col="#FF0000FF",signiffade="55",
                  xpd=NA,
                  text=c("Significant","Not Significant")
                  ){
    tpar=graphics::par(xpd=xpd);
    rcol=paste0(substr(col,1,7),signiffade)
    x=c(x,x)
    y=c(y,y+2*size[2]+gap)
    graphics::rect(x[1]-size[1],y[1]-size[2],x[1]+size[1],y[1]+size[2],
         col=imageback,border=border)
    graphics::rect(x[1]-size[1],y[1]-size[2],x[1]+size[1],y[1]+size[2],
         col=col,border=NA)
    graphics::rect(x[2]-size[1],y[2]-size[2],x[2]+size[1],y[2]+size[2],
         col=imageback,border=border)
    graphics::rect(x[2]-rectdelta[1]*size[1],y[2]-rectdelta[2]*size[2],
         x[2]+rectdelta[1]*size[1],y[2]+rectdelta[2]*size[2],
         col=rcol,border=NA)
    graphics::text(x[1]+size[1]+laboffset,y[1],text[1],adj=0,cex=cex.text)
    graphics::text(x[2]+size[1]+laboffset,y[2],text[2],adj=0,cex=cex.text)
    graphics::par(tpar)
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
#' @return an invisble list of the used plot symbols for each curve
#' @seealso \code{\link{plot.ClarityScan}} for the main plotting interface for persistences.
#' @examples
#' \donttest{
#' scan=Clarity_Scan(dataraw) ## Core Clarity
#' predmix=Clarity_Predict(datamix,scan) ## Core prediction
#' predrep=Clarity_Predict(datarep,scan) 
#' 
#' scanmix=Clarity_Scan(datamix)
#' scanrep=Clarity_Scan(datarep)
#' ## Plot the objectives
#' Clarity_ObjectivePlot(list(scan=scan,
#'                            scanmix=scanmix,
#'                            scanrep=scanrep,
#'                            predrep=predrep,
#'                            predmix=predmix),
#'                       col=c(1,2,3,2,3),
#'                       lty=c(1,1,1,2,2),
#'                       ylim=c(1e1,5e4),
#'                       las=1,
#'                       log="y")
#' }
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
                                                    col=tcol,lwd=tlwd,lty=tlty,pch=tpch,extern=FALSE,...)
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
