% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Clarity.R
\name{Clarity_Predict}
\alias{Clarity_Predict}
\title{Predict one Similarity from a Clarity or ClarityScan object}
\usage{
Clarity_Predict(Ynew, clist, Ysvd = NULL, ...)
}
\arguments{
\item{Ynew}{The data to be predicted}

\item{clist}{The learned Clarity object, either of class "Clarity" or of class "ClarityScan"}

\item{Ysvd}{Default NULL. SVD of Y, which is needed for the method="SVD". You can provide it or it will be computed for you when using this method. NB: It should be the SVD of Y-rowMeans(Y) to remove the mean effect, as computed by \code{\link{c_svd}}.}

\item{...}{Additional arguments to \code{\link{c_svdlist}}, if the original method used this.}
}
\value{
An object of the same class as clist provided, with updated Y, X and derived features.
}
\description{
Takes each Clarity object and learns a new X for the new data Ynew provided.
}
\examples{
\dontrun{
scanraw=Clarity_Scan(dataraw) # Generate an initial run
# Apply it to a new dataset with the same structure
scanrepfromraw=Clarity_Predict(datarep,clist=scanraw) 
# Apply it to a new dataset with slightly different structure
scanmixfromraw=Clarity_Predict(datamix,clist=scanraw) 
}
}
\seealso{
\code{\link{Clarity_Scan}} to generate an estimate for clist to be used here, \code{\link{plot.Clarity}} and \code{\link{plot.ClarityScan}} for plotting.
}
\keyword{mixture}
