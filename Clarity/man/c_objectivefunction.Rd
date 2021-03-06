% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Clarity.R
\name{c_objectivefunction}
\alias{c_objectivefunction}
\title{Default objective function for Clarity}
\usage{
c_objectivefunction(A, X, Y)
}
\arguments{
\item{A}{An N by K non-negative matrix of mixture proportions}

\item{X}{A K by K non-negative matrix of latent weights}

\item{Y}{An N by N non-negative matrix of similarities}
}
\value{
A numeric value of the objective function
}
\description{
This function evaluates
     || Y - A X A^T ||_2
}
\keyword{mixture}
