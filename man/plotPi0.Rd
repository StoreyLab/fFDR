\name{plotPi0}
\alias{plotPi0}
\title{Plot estimated pi0 as a function of qZ}
\usage{
plotPi0(fq, horizontal.line = TRUE, ...)
}
\arguments{
  \item{fq}{The data.table returned from fqvalue}

  \item{horizontal.line}{Whether to include a dashed
  horizontal line representing the traditional estimate of
  pi0}

  \item{...}{Extra arguments to give to qvalue for
  estimating traditional pi0}
}
\value{
ggplot object (print to show the plot)
}
\description{
Plot predicted pi0 as a function of the quantile of Z, and
optionally include a horizontal line for the traditional
pi0 estimate.
}

