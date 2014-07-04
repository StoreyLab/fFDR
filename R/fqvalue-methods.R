# methods for an fqvalue object (besides plots, which are handled in plots.R)

setMethod("length", "fqvalue", function(x) nrow(x@table))
setMethod("as.numeric", "fqvalue", function(x, ...) x@table$fqvalue)
setMethod("as.data.frame", "fqvalue", function(x, ...) as.data.frame(x@table))

setMethod("show", "fqvalue",
          function(object) {
              cat("Estimation of functional qvalue on", length(object), "pvalues\n")
              cat("Use plot() on this object to construct a z vs pvalue scatterplot.",
                  "Use as.numeric() on this object to access the vector of fqvalues,",
                  "or as.data.frame() to extract a table comparing",
                  "p-values, z, and fqvalue\n")
          })



## Plots

#' Scatterplot comparing the p-value to the quantile of Z
#' 
#' Construct a scatterplot comparing the p-value to the quantile of Z,
#' coloring points either based on an oracle TRUE/FALSE (if the column
#' \code{oracle} is included) or based on thresholds of the functional
#' qvalue.
#' 
#' @param x fqvalue object
#' @param threshold vector of thresholds used for colors on the graph
#' (default c(.005, .01, .05, .1))
#' @param cloglog whether the plot should be on a log(-log(pvalue)) scale
#' (default FALSE)
#' @param pi0 Should a second panel of the estimate of pi0(z) be shown
#' underneath the scatterplot
#' @return ggplot object (print to show the plot)
#' 
#' @import ggplot2
#' @import scales
#' @import reshape2
#' @import data.table
#' @import gridExtra
#' 
#' @export
setMethod("plot", c("fqvalue"),
          function(x, pi0=TRUE, threshold=c(.005, .01, .05, .1),
                   cloglog=FALSE, ...) {
              tab = copy(x@table)
              qv = qvalue(tab$pvalue)$qvalue
              
              num.below = sapply(threshold, function(q) sum(qv < q))
              vlines = data.frame(qvalue=as.factor(threshold),
                                  pvalue=c(0, sort(tab$pvalue))[num.below + 1])
              qp = factor(c(threshold, 1))
              tab$confidence = qp[findInterval(tab$fqvalue, threshold) + 1]

              if (!any(tab$fqvalue < max(threshold))) {
                  # if nothing is significant, just show everything
                  highest.pval = 1
              } else {
                  highest.pval = max(tab$pvalue[tab$fqvalue < max(threshold)])
              }
              tab = tab[tab$pvalue <= max(highest.pval, max(vlines$pvalue))]
              
              g = (ggplot(tab, aes(x=z, y=pvalue)) +
                       scale_color_brewer(palette="Spectral") +
                       theme(axis.text.x = element_text(angle = 30, hjust = 1)))
              
              if (cloglog) {
                  g = g + scale_y_continuous(trans=reverseloglog_trans(10))
                  pval.points = 10^(-10^(seq(-4.5, 2, .05)))
              }
              else {
                  pval.points = seq(0.002, 1, .002)
              }
              
              qZ.points = seq(0, 1, .001)
              
              if (!is.null(tab$oracle)) {
                  g = g + geom_point(aes(col=oracle), data=tab, size=3)
              }
              else {
                  g = g + geom_point(aes(col=confidence),
                                     data=tab, size=3)
              }
              # , color=factor(..level..))
              #g = g + stat_contour(aes(x=X2, y=X1, z=value), data=qval.m, breaks=threshold)
              if (is.null(tab$oracle)) {
                  g = g + geom_hline(aes(yintercept=pvalue, col=qvalue), data=vlines)
              }
              #scale_color_gradient(name="Q-value", low="red", high="white") +

              if (pi0) {
                  g = arrangeGrob(g, plot(x@fPi0), nrow=2)
              }

              return(g)
          })



#' Compare functional q-value to traditional q-values
#' 
#' Compare traditional qvalue to functional q-values in a scatterplot,
#' colored either by the quantile of Z or by a given oracle
#' 
#' @param fq An fqvalue object
#' @param ... additional arguments passed to the traditional qvalue function
#' 
#' @return ggplot object (print to show the plot)
#' 
#' @details If the column \code{oracle} exists in the fqvalue, compareQvalue
#' will use that as the colors of the points. Otherwise, the color depends
#' on the quantiles of Z.
#' 
#' @import ggplot2
#' @import qvalue
#' 
#' @export
compareQvalue = function(fq, ...) {
    tab = fq@table
    tab$qvalue = qvalue(tab$pvalue, ...)$qvalue
    g = (ggplot(tab, aes(qvalue, fqvalue,
                         label=round(pvalue, 3))) + xlab("Traditional q-value") +
             ylab("Functional q-value"))
    if (is.null(tab$oracle)) {
        g = g + geom_point(aes(col=z))
    }
    else {
        g = g + geom_point(aes(col=oracle))
    }
    g + geom_abline(col="red")
}


#' Plot curves of local FDR or density (f) by p-value
#' 
#' Plot predicted local FDR as a function of p-value, separate curves for
#' different qZ values. Shown 
#' 
#' @param fq The data.table returned from fqvalue
#' @param cutoff fqvalue cutoff at which to draw horizontal line
#' 
#' @return ggplot object (print to show the plot)
#' 
#' @import reshape2
densityCurves = function(fq, cutoff=.05, plottype="lfdr", ...) {
    stop("Currently Deprecated")
    library(ggplot2)
    library(qvalue)
    
    plottype = match.arg(plottype, c("lfdr", "density"))
    
    density.grid = attr(fq, plottype)
    m = as.data.table(melt(density.grid))
    m$qZ = k$x[m$Var1]
    m$pvalue = k$pvalue[m$Var2]
    
    g = (ggplot(m[Var2 %% 5 == 0 & Var1 %% 20 == 0], aes(pvalue, value, col=qZ, group=qZ)) +
             geom_line())
    g
}


#' log(-log(x)) transformation for a ggplot2 scale
#' 
#' @description this code is a modified version of this:
#' http://stackoverflow.com/questions/11053899
reverseloglog_trans <- function(base = exp(1)) {
    trans <- function(x) log(-log(x, base), base)
    inv <- function(x) base^(-(base^x))
    trans_new(paste0("reverseloglog-", format(base)), trans, inv, 
              function(x) c(.9985, .99, .93, .65, .1, 1e-05,
                            1e-25, 1e-100), 
              domain = c(1e-300, Inf))
}

#' Plot demonstrating how lambda is chosen to minimize MISE in
#' estimation of pi0(z)
#' 
#' Performs a plotMISE plot of the functional pi0 estimate contained
#' within an fqvalue object
#' 
#' @param x fqvalue object
#' 
#' @export
setMethod("plotMISE", "fqvalue", function(x, ...) plotMISE(x@fPi0))
    