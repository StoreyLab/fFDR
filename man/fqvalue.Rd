% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fqvalue.R
\name{fqvalue}
\alias{fqvalue}
\title{Estimate functional q-value based on p-value and an informative variable}
\usage{
fqvalue(
  p.value,
  z0,
  pi0.method = "gam",
  lambda = seq(0.4, 0.9, 0.1),
  fixed.pi0 = FALSE,
  monotone.window = 0.01,
  ...
)
}
\arguments{
\item{p.value}{A vector of p-values}

\item{z0}{A vector of observations from the informative variable, of the same length as \code{p}}

\item{pi0.method}{Method for estimating the functional proportion pi0(z); either "gam" (default), "glm", "kernel" or "bin"}

\item{lambda}{The parameter Lambda or its values to use in estimating pi0(z)}

\item{fixed.pi0}{Whether pi0(z) is believed to be independent of z}

\item{monotone.window}{Parameter used to force estimated densities to
decrease with increasing p-values; higher value of this parameter means densities are smoothed
more aggressively. If NULL, perform no such smoothing}

\item{...}{Extra arguments to be passed to kernelUnitInterval for estimating the density}
}
\value{
An object of S3 class "fqvalue"
}
\description{
Estimate functional q-values based on p-values and realizations of
the informative variable z, where z may affect either the power of a statistical test or the likelihood of a true null hypothesis.
}
\details{
Assume the random variable z0 may affect the power of a statistical test (that induces the p-values) or the likelihood of a true null hypothesis. The m observations z0_i, i=1,...,m of z0 are quantile transformed into z_i, i=1,...,m such that z_i = rank(z0_i) / m, where rank(z0_i) is the rank of z0_i among z0_i, i=1,...,m. Consequently, z_i, i=1,...,m are approximately uniformly distributed on the interval [0,1]. When z_i, i=1,...,m are regarded as observations from the random variable z, then z is approximately uniformly distributed on [0,1]. Namely, z0 has been quantile transformed into z, and they are equivalent. Further, z or z0 is referred to as the informative variable.

The likelihood of a true null hypothesis is regarded as a function of z, referrred to as the functional proportion, and denoted by pi0(z), and the fFDR methodology is applied to the m paired observations (p_i,z_i), i=1,...,m of the p-value p and z.
}
