\name{fqvalue}
\alias{fqvalue}
\title{Estimate functional q-values based on p-values and a power surrogate}
\usage{
fqvalue(p, Z, oracle = NULL, fixed.pi0 = FALSE, p.bandwidth = NULL,
  z.bandwidth = NULL, grid.points = 250, transformation = "cloglog",
  min.pval = 1e-30, only.pi0 = FALSE, ...)
}
\arguments{
  \item{p}{A vector of p-values}

  \item{Z}{a vector of the power surrogate variable, where
  a higher Z generally represents a higher power of the
  hypothesis test.}

  \item{fixed.pi0}{Whether the true pi0 (null hypothesis
  rate) is believed to be the same for all Z.}

  \item{lambda}{The value of the tuning parameter to
  estimate pi0.}

  \item{pi0.method}{Method used to estimate pi0 as a
  function of Z: either "kernel" (default), "spline",
  "logistic", or "LOESS". Ignored if fixed.pi0 is TRUE}

  \item{bandwidth}{The scaling bandwidth to use for the
  kernel density estimation: either a single value or a
  vector of bandwidths for smoothing in the qZ and p-value
  directions, respectively.}

  \item{grid.points}{The number of grid points to use in
  each direction in kernel smoothing, either a single value
  or a vector of integers for smoothing in the qZ and
  p-value directions, respectively.}

  \item{transformation}{Transformation of p-values before
  the kernel density is estimated: either "cloglog"
  (default)=log(-log(p)) or "probit", the inverse of the
  cumulative density of the normal distribution.}

  \item{pi0.bandwidth}{If pi0.method is "kernel" (ignored
  otherwise), the scaling bandwidth to use for the kernel
  estimation of pi0}

  \item{smooth.df}{If pi0.method is "spline" (ignored
  otherwise), the number of degrees of freedom to use. If
  NULL (default), use cross validation to decide}

  \item{only.pi0}{Return as soon as (f)pi0 is calculated,
  without calculating fqvalue. Useful in simulation}

  \item{...}{Extra arguments to be passed to either}
}
\value{
A data.table with columns:

\item{pvalue}{The original pvalues given to the function
(in the same order)} \item{Z}{The original Z given to the
function} \item{qZ}{The quantiles of Z, which are used for
the nonparametric calculations of pi0 and fqvalue.}
\item{pi0}{Functional pi0, which depends on qZ.}
\item{fqvalue}{Functional q-value: the estimated false
discovery rate of rejecting all hypotheses with this
f-qvalue or lower}
}
\description{
Estimate functional q-values based on a vector of p-values
and a "power surrogate" Z, where Z is a known value that
tends to increase with the power of the test.
}

