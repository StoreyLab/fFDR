\name{scatterZPvalue}
\alias{scatterZPvalue}
\title{Scatterplot comparing the p-value to the quantile of Z}
\usage{
scatterZPvalue(fq, threshold = c(0.005, 0.01, 0.05, 0.1), doublelog = FALSE,
  pi0 = FALSE)
}
\arguments{
  \item{fq}{fqvalue object}

  \item{threshold}{vector of thresholds used for colors on
  the graph (default c(.005, .01, .05, .1))}

  \item{doublelog}{whether the plot should be on a
  log(-log(pvalue)) scale (default FALSE)}

  \item{pi0}{Should a line plot of pi0 be shown underneath
  the plot}
}
\value{
ggplot object (print to show the plot)
}
\description{
Construct a scatterplot comparing the p-value to the
quantile of Z, coloring points either based on an oracle
TRUE/FALSE (if the column \code{oracle} is included) or
based on thresholds of the functional qvalue.
}

