\name{getZifNetwork}
\alias{getZifNetwork}
\title{Put coefficients from a set of network regressions into a matrix/array}
\usage{
  getZifNetwork(listOfFits, constraint = 1,
    collapse = FALSE, type = "grp.l1", norm.grid)
}
\arguments{
  \item{listOfFits}{output from fitZifNetwork}

  \item{constraint}{a constraint, which will be translated
  into the regularization tuning parameter lambda for each
  regression}

  \item{collapse}{should a collapsed, symmetrized, boolean
  adjacency matrix be returned?}

  \item{type}{the type of constraint imposed}

  \item{norm.grid}{optional output from fortify.zifnetwork
  to speed calculation}
}
\value{
  an array
}
\description{
  The array is ngenes X ngenes X {2,4}, with the last
  dimension depending on whether zero-inflated or hurdle
  predictors were used. When type = 'grp.l1' lambda for
  each regression is select to minimize the hurdle deviance
  subject to the grouped l1 norm being less than
  'constraint'. When type = 'comb.l1', ditto but with the
  regular l1 norm being constrained. When type =
  'lambda.min' lambda is selected to be 'constraint' times
  the lambda.min determined by cross validation.
}

