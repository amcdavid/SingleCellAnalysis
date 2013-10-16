\name{plot.threshold}
\alias{plot.threshold}
\title{Plot NanoString log Counts vs cluster ID}
\usage{
  plot.threshold(M, densitiesFromThres, primerids)
}
\arguments{
  \item{M}{output from thresholdNanostring, debug=TRUE}

  \item{densitiesFromThres}{output from
  thresholdNanostring, debug=TRUE}

  \item{primerids}{character vector of primerids to plot}
}
\value{
  ggplot object
}
\description{
  Histogram of log Counts by primerid, colored by cluster
}
\author{
  andrew
}

