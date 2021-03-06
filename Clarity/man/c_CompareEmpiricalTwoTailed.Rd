% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ClarityBootstrap.R
\name{c_CompareEmpiricalTwoTailed}
\alias{c_CompareEmpiricalTwoTailed}
\title{Compare Observed to Bootstrap values (Empirical two tailed test)}
\usage{
c_CompareEmpiricalTwoTailed(x, o)
}
\arguments{
\item{x}{the distribution of theoretical or simulated values}

\item{o}{the observed value}
}
\value{
a p value
}
\description{
Compare an observation (o) with a set of values (x) and report the empirical p value for a one tailed test that o is far from the mean of x.

If the number of observations is low, consider \code{\link{c_CompareNormalTwoTailed}} instead.

For a one-tailed test, consider \code{\link{c_CompareEmpirical}} instead.
}
\seealso{
\code{\link{c_pval}} in which this function can be used.
}
