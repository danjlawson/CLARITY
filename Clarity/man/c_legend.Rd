% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ClarityPlots.R
\name{c_legend}
\alias{c_legend}
\title{Plot a legend that describes how CLARITY shows significance}
\usage{
c_legend(
  x,
  y,
  size = c(0.5, 0.5),
  rectdelta = c(0.8, 0.8),
  imageback = "white",
  border = "grey",
  cex.text = 1.5,
  gap = 0,
  laboffset = 0.1,
  col = "#FF0000FF",
  signiffade = "55",
  xpd = NA,
  text = c("Significant", "Not Significant")
)
}
\arguments{
\item{x}{the x location of the legend (center of the top element); you can specifiy top and bottom elements locations by giving two values.}

\item{y}{the y location of the legend (center of the top element); you can specifiy top and bottom elements locations by giving two values.}

\item{size}{(default=c(0.5,0.5)) the size of the element as a distance from the center specified as (x, y)}

\item{rectdelta}{(default=c(0.8,0.8)) the relative reduction of the non-significant rectangle}

\item{imageback}{(default="white") colour to show behind the non-significant rectangle}

\item{border}{(default="grey") colour of the border to place around each rectangle. Set to NA to omit.}

\item{cex.text}{(default=1.5) size of the legend text}

\item{gap}{(default=0) gap between the two elements, if their locations were not both specified}

\item{laboffset}{(default=0.1) distance between the edge of the rectangle and the start of the text}

\item{col}{(default="#FF0000FF") base colour for significance}

\item{signiffade}{(default="55") change to the alpha channel for non-significance}

\item{xpd}{(default=NA) how to protect plotting outside of the region; see \code{\link{par}}(xpd).}

\item{text}{(default=c("Significant","Not Significant")) what to write as the legends.}
}
\description{
This is a rather ugly, manual hack to add a legend outside of the main image area. It puts example rectangles that are significant and non-significant at a given location in the image co-ordinates.

You might want to show this some other way, but this can work with manual tuning of x,y, and size in image units.
}
\examples{
\dontrun{
scan=Clarity_Scan(dataraw) ## Core Clarity

## Plot persistence chart with p-values
scanpred=Clarity_Compare(scan,D=datarawD,Dnew=datamixD)
predplot=plot(scanpred)
## Add the legend outside of the regular image space
c_legend(0,22,size=c(10,0.5))
}
}
