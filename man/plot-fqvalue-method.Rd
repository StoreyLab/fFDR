\docType{methods}
\name{plot,fqvalue-method}
\alias{plot,fqvalue-method}
\title{Scatterplot comparing the p-value to the quantile of Z}
\usage{
\S4method{plot}{fqvalue}(x, pi0 = TRUE, threshold = c(0.005, 0.01, 0.05,
  0.1), cloglog = FALSE, ...)
}
\arguments{
  \item{x}{fqvalue object}

  \item{threshold}{vector of thresholds used for colors on
  the graph (default c(.005, .01, .05, .1))}

  \item{cloglog}{whether the plot should be on a
  log(-log(pvalue)) scale (default FALSE)}

  \item{pi0}{Should a second panel of the estimate of
  pi0(z) be shown underneath the scatterplot}

  \item{...}{Arguments to be passed to methods, such as
  \link{graphical parameters} (see \code{\link{par}}).
  Many methods will accept the following arguments:
  \describe{ \item{\code{type}}{what type of plot should be
  drawn.  Possible types are \itemize{ \item \code{"p"} for
  \bold{p}oints, \item \code{"l"} for \bold{l}ines, \item
  \code{"b"} for \bold{b}oth, \item \code{"c"} for the
  lines part alone of \code{"b"}, \item \code{"o"} for both
  \sQuote{\bold{o}verplotted}, \item \code{"h"} for
  \sQuote{\bold{h}istogram} like (or \sQuote{high-density})
  vertical lines, \item \code{"s"} for stair \bold{s}teps,
  \item \code{"S"} for other \bold{s}teps, see
  \sQuote{Details} below, \item \code{"n"} for no plotting.
  } All other \code{type}s give a warning or an error;
  using, e.g., \code{type = "punkte"} being equivalent to
  \code{type = "p"} for S compatibility.  Note that some
  methods, e.g. \code{\link{plot.factor}}, do not accept
  this.  } \item{\code{main}}{an overall title for the
  plot: see \code{\link{title}}.} \item{\code{sub}}{a sub
  title for the plot: see \code{\link{title}}.}
  \item{\code{xlab}}{a title for the x axis: see
  \code{\link{title}}.} \item{\code{ylab}}{a title for the
  y axis: see \code{\link{title}}.} \item{\code{asp}}{the
  \eqn{y/x} aspect ratio, see \code{\link{plot.window}}.} }
  }
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

