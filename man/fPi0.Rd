\name{fPi0}
\alias{fPi0}
\title{estimate pi0 as a function of Z}
\usage{
fPi0(p, Z, lambda = seq(0.3, 0.9, 0.1), method = "glm", df = 3,
  breaks = 5, transformation = "cloglog", bs.its = 100, ...)
}
\arguments{
  \item{p}{a vector of p-values}

  \item{Z}{a factor on which pi0 is believed to depend
  \code{p}}

  \item{lambda}{}

  \item{method}{Either "glm" (default), "gam", "kernel", or
  "bin"}

  \item{df}{Degrees of freedom to use in spline in "gam"
  method}

  \item{transformation}{Transformation to use in "kernel"
  method, either "cloglog" (default), "probit", or "ident"}

  \item{breaks}{Either a number of (evenly spaced) break
  points in "bin" method, or a vector of break points (from
  0 to 1) to use for bins}

  \item{...}{Additional arguments to glm, gam or the kernel
  estimator}
}
\value{
a vector of estimated pi0
}
\description{
\code{fpi0} estimates a pi0 function that depends on
another parameter Z
}
\details{
In short, the "glm", "gam", and "kernel" methods attempt to
estimate:

\code{pi0(z) = Pr(p>lambda|z)/(1-lambda)}

The glm and gam approaches define a variable
\code{eta=I{p>lambda}}, and use a modification of logistic
regression to fit \code{eta~f(z)}. The kernel density
estimate examines the density of z where p>lambda and
computes Pr(p>lambda) from that.

Binning simply computes the Storey estimate with the given
lambda within each bin.
}

