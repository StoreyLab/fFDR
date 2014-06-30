\name{fixedPi0}
\alias{fixedPi0}
\title{compute pi0 where the truth is independent of Z}
\usage{
fixedPi0(p, Z, lambda = seq(0, 0.9, 0.05), tau = c(0, 0.9, 0.05))
}
\arguments{
  \item{p}{a vector of p-values}

  \item{Z}{a vector of a power-informative factor, of the
  same length as \code{p}}

  \item{lambda}{Possible values of lambda to try}

  \item{tau}{Possible values of tau to try}
}
\value{
estimated pi0
}
\description{
\code{fixedPi0} estimates pi0 for a set of p-values using a
power-informative factor Z. This the assumption that
although power depends on Z, the truth of the hypotheses is
independent of it.
}

