\name{kernelUnitInterval}
\alias{kernelUnitInterval}
\title{estimate density on the unit interval; either 1d or 2d, using
local regression desity estimation}
\usage{
kernelUnitInterval(x, transformation = "probit", eval.points = x,
  cv = TRUE, ...)
}
\arguments{
  \item{x}{either a vector or a 2-column matrix}

  \item{transformation}{either probit (default),
  complementary log-log, or identity (not recommended)}

  \item{eval.points}{Points at which to evaluate the
  estimate, default x}

  \item{cv}{Whether to use generalized cross-validation to
  choose the nn (nearest neighbor) smoothing parameter}

  \item{...}{additional arguments to be passed to lp in
  locfit, used only if cv=FALSE}
}
\description{
Used in the fpi0 and fqvalue density estimation steps
}

