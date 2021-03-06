% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimation.R
\name{mle_cens_gpd}
\alias{mle_cens_gpd}
\title{A function that calculates the maximum likelihood estimator based on the censored generalized Pareto
distribution.}
\usage{
mle_cens_gpd(data)
}
\arguments{
\item{data}{Dataframe. One column is called "Injury_Length" and corresponds to the observations. One
column called "Censored" of 0s and 1s. 0 indicates that the observation is not censored. 1 indicates the
observation is censored.}
}
\value{
Numeric. Corresponding to the maximum likelihood estimator based on the censored generalized
Pareto distribution.
}
\description{
A function that calculates the maximum likelihood estimator based on the censored generalized Pareto
distribution.
}
\examples{

data <- data.frame(Injury_Length = rexp(10), Censored = rbinom(10, 1, 0.5))

mle_cens_gpd(data = data)
}
\author{
Annabelle Redelmeier
}
