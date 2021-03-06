% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ClaritySim.R
\name{claritysim_dist}
\alias{claritysim_dist}
\title{Standard Distance function for Clarity}
\usage{
claritysim_dist(x)
}
\arguments{
\item{x}{An N by K matrix of N subjects observed at K features}
}
\value{
An N by N matrix of the distances
}
\description{
Computes the euclidean distance between all points using as.matrix(dist(x)), enforcing the diagonal to take the minimum value of the rest of each row.
}
\seealso{
\code{\link{resimulatedDistances}}, which computes this "correctly" in a way that is unbiased for the diagonal.
}
