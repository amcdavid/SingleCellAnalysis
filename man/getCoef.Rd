% Generated by roxygen2 (4.0.1): do not edit by hand
\name{getCoef}
\alias{getCoef}
\title{Get coefficients and standard errors from a zlm fit}
\usage{
getCoef(zlm, cname)
}
\arguments{
\item{zlm}{zlm object}

\item{cname}{character vector of coefficient names}
}
\value{
3D array: (disc, cont) X (coefficient, standard error) X coefficient
}
\description{
For each coefficient name in cname, return discrete and continuous values and standard errors
}

