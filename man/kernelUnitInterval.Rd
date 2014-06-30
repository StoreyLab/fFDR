\name{kernelUnitInterval}
\alias{kernelUnitInterval}
\title{estimate density on the unit interval; either 1d or 2d, using
local regression desity estimation}
\usage{
kernelUnitInterval(x, transformation = "probit", eval.points = x,
  cv = TRUE, epsilon = 1e-15, maxk = 100, grid = FALSE, ...)
}
\arguments{
  \item{x}{either a vector or a 2-column matrix}

  \item{transformation}{either probit (default),
  complementary log-log, or identity (not recommended)}

  \item{eval.points}{Points at which to evaluate the
  estimate, default x}

  \item{cv}{Whether to use generalized cross-validation to
  choose the nn (nearest neighbor) smoothing parameter.
  Currently not permitted for 2d estimation}

  \item{epsilon}{How close are allowed to come to 0 or 1;
  values are thresholded to be in the range [epsilon,
  1-epsilon]}

  \item{maxk}{maxk argument passed to locfit}

  \item{...}{additional arguments to be passed to lp in
  locfit, used only if cv=FALSE}
}
\description{
Used in the fpi0 and fqvalue density estimation steps
}

