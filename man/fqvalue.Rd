\name{fqvalue}
\alias{fqvalue}
\title{Estimate functional q-values based on p-values and an informative factor
z0}
\usage{
fqvalue(pvalue, z0, fixed.pi0 = FALSE, monotone.window = 0.01,
  transformation = "cloglog", ...)
}
\arguments{
  \item{pvalue}{A vector of p-values}

  \item{z0}{a vector providing information either on the
  power of each test or the likelihood of the null}

  \item{fixed.pi0}{Whether the true pi0 (null hypothesis
  rate) is believed to be the same for all z.}

  \item{monotone.window}{Parameter used to force estimated
  densities to decrease with increasing p-values- higher
  means densities are smoothed more aggressively. If NULL,
  perform no such smoothing.}

  \item{transformation}{Transformation of p-values before
  the kernel density is estimated: either "cloglog"
  (default)=log(-log(p)) or "probit", the inverse of the
  cumulative density of the normal distribution.}

  \item{...}{Extra arguments to be passed to estFPi0}
}
\value{
An object of class \linkS4class{fqvalue}
}
\description{
Estimate functional q-values based on a vector of p-values
and z0, where z0 is known to be related to either power or
the likelihood of the null hypothesis.
}

