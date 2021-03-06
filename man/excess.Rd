% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper_functions.R
\name{excess}
\alias{excess}
\title{Calculate excess above the threshold. Threshold calculated using ne = count of excesses above the threshold.}
\usage{
excess(data, ne = NULL, threshold = NULL)
}
\arguments{
\item{data}{Dataframe. Contains at least the columns "ID" and "Injury_Length"..}

\item{ne}{Integer. Count of excesses above the threshold.}

\item{threshold}{Integer. If \code{ne} is not supplied, \code{threshold} can be used to specify the threshold
to use.}
}
\value{
Dataframe. Of observations above the threshold (decided by \code{ne}) and the excess above the threshold.
}
\description{
Calculate excess above the threshold. Threshold calculated using ne = count of excesses above the threshold.
}
\examples{
data <- data.frame(ID = c(rep(1, 3), rep(2, 2), rep(3, 4), 5), Injury_Length = rexp(10))

excess(data, ne = 1)
}
\author{
Annabelle Redelmeier
}
