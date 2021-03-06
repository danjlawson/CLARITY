% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ClarityMultiplicative.R
\name{c_DecreaseK}
\alias{c_DecreaseK}
\title{Make a smaller Clarity object by randomly peturbing a larger one}
\usage{
c_DecreaseK(Y, A, X)
}
\arguments{
\item{Y}{Non-negative Similiarty matrix being fitted (N by N)}

\item{A}{Non-negative mixture matrix}

\item{X}{Non-negative latent similarity matrix}
}
\value{
A list containing the following:
\itemize{
\item X A K by K non-negative matrix of latent weights
\item A An N by K non-negative matrix of mixture proportions
\item dist The objective function evaluated at the new A and X
}
}
\description{
Make a Clarity object of size K-1 from one of size K-1 by removing a random column from of size K
}
\keyword{mixture}
