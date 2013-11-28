\name{pearson.pairwise}
\alias{pearson.pairwise}
\title{Calculate pairwise correlation coefficient}
\usage{
  pearson.pairwise(sc, diagonals = "zero", partial = FALSE)
}
\arguments{
  \item{sc}{SingleCellAssay object}

  \item{diagonals}{currently ignored}

  \item{partial}{adjust for ngeneson}
}
\value{
  a PairwisePearson object, which is a matrix with slots
  "nsamp" giving the number of cells the cor. coef. were
  calculated over, "npos" giving the number of expressed
  cells in each gene.
}
\description{
  Calculates the correlation coefficient between each pair
  of genes on the dichotomous level. The diagonal is set to
  zero to avoid having the fisher transform blow up.
}

