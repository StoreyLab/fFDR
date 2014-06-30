\name{densityCurves}
\alias{densityCurves}
\title{Plot curves of local FDR or density (f) by p-value}
\usage{
densityCurves(fq, cutoff = 0.05, plottype = "lfdr", ...)
}
\arguments{
  \item{fq}{The data.table returned from fqvalue}

  \item{cutoff}{fqvalue cutoff at which to draw horizontal
  line}
}
\value{
ggplot object (print to show the plot)
}
\description{
Plot predicted local FDR as a function of p-value, separate
curves for different qZ values. Shown
}

