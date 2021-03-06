% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ClaritySVD.R
\name{c_SVD_ResidualExtendFixedK}
\alias{c_SVD_ResidualExtendFixedK}
\title{Extend a Clarity result by fitting to the residuals, using the SVD method and a fixed k}
\usage{
c_SVD_ResidualExtendFixedK(clist, Rsvd, R, kextra, verbose = TRUE)
}
\arguments{
\item{clist}{the Clarity object that was learned, as return by \code{\link{Clarity_Predict}}}

\item{Rsvd}{The SVD of R-rowMeans(R).}

\item{R}{The distance matrix of the residuals, to be fit}

\item{kextra}{The number of additional mixture components}

\item{verbose}{(default TRUE) Whether to output progress information}
}
\value{
A ClarityExtend object as described in \code{\link{Clarity_Extend}}
}
\description{
Given a Clarity object, and details of the residuals that it has failed to fit, fit an additional kextra mixture components to it to learn its unique model structure.  This uses the spectral embedding of the original fit, and the spectral embedding of the residuals, to make a feature vector that can be separated into the two components.
}
\seealso{
\code{\link{Clarity_Extend}} for the recommended interface.
}
