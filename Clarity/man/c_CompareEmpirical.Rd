% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ClarityBootstrap.R
\name{c_CompareEmpirical}
\alias{c_CompareEmpirical}
\title{Compare Observed to Bootstrap values (Empirical one tailed test)}
\usage{
c_CompareEmpirical(x, o)
}
\arguments{
\item{x}{the distribution of theoretical or simulated values}

\item{o}{the observed value}
}
\value{
a p value
}
\description{
Compare an observation (o) with a set of values (x) and report the empirical p value for a one tailed test that o is large compared to x.

If the number of observations is low, consider \code{\link{c_CompareNormal}} instead.

For a two-tailed test, consider \code{\link{c_CompareEmpiricalTwoTailed}} instead.
}
\seealso{
\code{\link{c_pval}} in which this function can be used.
}
