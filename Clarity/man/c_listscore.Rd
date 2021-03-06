% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Clarity.R
\name{c_listscore}
\alias{c_listscore}
\title{Score a list of Clarity objects}
\usage{
c_listscore(clist)
}
\arguments{
\item{clist}{A list of Clarity objects as returned by \code{\link{Clarity_Scan}}}
}
\value{
A numeric vector of length equal to the length of clist, containing the objective function evaluations
}
\description{
Takes a list of Clarity results and evaluate their objective function at the current settings
}
\keyword{mixture}
