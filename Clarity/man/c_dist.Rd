% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ClarityCompare.R
\name{c_dist}
\alias{c_dist}
\title{Standard Distance function for Clarity}
\usage{
c_dist(x)
}
\arguments{
\item{x}{An N by K matrix of N subjects observed at K features}
}
\value{
An N by N matrix of the distances
}
\description{
Computes the euclidean distance between all points using as.matrix(dist(x)), using \code{\link{c_fixdiagonals}} to enforce the diagonal to take the minimum value of the rest of each row.

Properly, the diagonal should be be the distance that would be achieved by a replicate of the data x, but this is typically unavailable without making assumptions on the distribution of the features. Using the minimum is typically appropriate; see simulationToDistance in the package ClaritySim for details and the help from that package for a demonstration.
}
\seealso{
\code{\link{Clarity_Compare}}, which takes an argument distfn which can be given this or another function that generates a similarity from feature data.
}
