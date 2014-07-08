\name{fqvalue}
\alias{fqvalue}
\title{Estimate functional q-values based on p-values and an informative factor
z0}
\usage{
fqvalue(pvalue, z0, pi0.method = "gam", lambda = seq(0.4, 0.9, 0.1),
  fixed.pi0 = FALSE, monotone.window = 0.01, transformation = "probit",
  ...)
}
\arguments{
  \item{pvalue}{A vector of p-values}

  \item{z0}{a vector providing information either on the
  power of each test or the likelihood of the null}

  \item{pi0.method}{Method for estimating pi0(z): either
  "gam" (default)}

  \item{lambda}{Lambda parameter, or choices of parameter,
  to use in estimating pi0(z)}

  \item{fixed.pi0}{Whether the true pi0 (null hypothesis
  rate) is believed to be the same for all z.}

  \item{monotone.window}{Parameter used to force estimated
  densities to decrease with increasing p-values- higher
  means densities are smoothed more aggressively. If NULL,
  perform no such smoothing.}

  \item{...}{Extra arguments to be passed to
  kernelUnitInterval for estimating the density}
}
\value{
An object of class \linkS4class{fqvalue}
}
\description{
Estimate functional q-values based on a vector of p-values
and z0, where z0 is known to be related to either power or
the likelihood of the null hypothesis.
}

