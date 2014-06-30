# FPi0 plots

#' Given an fPi0 object, plot the estimates of pi_0(z) with multiple
#' choices of lambda
#' 
#' @param x fPi0 object
#' @param ... other plotting parameters (not used)
#' 
#' @return a ggplot2 graph plotting pi0 as a function of z
#' 
#' @import ggplot2
#' 
#' @export
setMethod("plot", c("fPi0"),
    function(x, ...) {
        x@tableLambda$Chosen = (x@tableLambda$lambda == x@lambda)
        g = ggplot(x@tableLambda, aes(z, fpi0, col=lambda, group=lambda,
                                  lty=Chosen)) + geom_line()
        g = g + scale_linetype_manual(values=c(3, 1))
        g = g + labs(col=expression(lambda))
        g = g + ylab(expression(hat(pi)[0]^(lambda) ~ (z)))
        return(g)
    })


#' Plot demonstrating how lambda is chosen to minimize MISE
#' 
#' Construct a plot containing two panels. The first compares the curve
#' estimated from each choice of lambda to a reference obtained from
#' the lowest lambda value (used to estimate how accurate the shape is).
#' The second panel shows how omega (error of the shape), delta.sq
#' (squared bias), and their sum (estimate of MISE) change with changing
#' lambda. The software's choice of lambda is highlighted as a vertical
#' dashed line.
#' 
#' @param x fpi0 object
#' 
#' @import ggplot2
#' @import reshape2
#' @import gridExtra
#' 
#' @export
setMethod("plotMISE", "fPi0", function(x, ...) {
    cm = melt(x@tableLambda[, !"k", with=FALSE], id=c("lambda", "z"))

    # plot comparing fpi0 to the reference at each stage
    cm$lambda = paste(cm$lambda, ifelse(cm$lambda == x@lambda, "(Chosen)", ""))
    cm$variable = ifelse(cm$variable == "fpi0", "Estimate", "Reference")
    g1 = ggplot(cm, aes(z, value, lty=variable)) + geom_line() +
        facet_wrap(~ lambda)
    g1 = g1 + labs(lty="")
    g1 = g1 + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    # plot comparing the bias, variance and MISE estimates
    sm = melt(x@MISE, id="lambda")
    g2 = ggplot(sm, aes(lambda, value, col=variable)) + geom_line()
    g2 = g2 + geom_vline(xintercept=x@lambda, lty=2)
    g2 = g2 + xlab(expression(lambda)) + labs(col="")
    
    print(arrangeGrob(g1, g2, nrow=2))
})


## fqvalue plots

#' Scatterplot comparing the p-value to the quantile of Z
#' 
#' Construct a scatterplot comparing the p-value to the quantile of Z,
#' coloring points either based on an oracle TRUE/FALSE (if the column
#' \code{oracle} is included) or based on thresholds of the functional
#' qvalue.
#' 
#' @param fq fqvalue object
#' @param threshold vector of thresholds used for colors on the graph
#' (default c(.005, .01, .05, .1))
#' @param doublelog whether the plot should be on a log(-log(pvalue)) scale
#' (default FALSE)
#' @param pi0 Should a line plot of pi0 be shown underneath the plot
 #' 
#' @return ggplot object (print to show the plot)
#' 
#' @export
scatterZPvalue = function(fq, threshold=c(.005, .01, .05, .1),
                        doublelog=FALSE, pi0=FALSE) {
    library(scales)
    library(reshape)
    library(data.table)
    
    qv = qvalue(fq$pvalue)$qvalue

    num.below = sapply(threshold, function(q) sum(qv < q))
    vlines = data.frame(qvalue=as.factor(threshold),
                        pvalue=sort(fq$pvalue)[num.below])
    qp = factor(c(threshold, 1))
    fq$confidence = qp[findInterval(fq$fqvalue, threshold) + 1]
    
    if (doublelog) {
        # have to get rid of the few highest values
        fq = fq[fq$pvalue < .999]
    }
    
    g = (ggplot(fq, aes(x=qZ, y=pvalue)) +
             scale_color_brewer(palette="Spectral") +
             theme(axis.text.x = element_text(angle = 30, hjust = 1)))
    
    if (doublelog) {
        g = g + scale_y_continuous(trans=reverseloglog_trans(10))
        pval.points = 10^(-10^(seq(-4.5, 2, .05)))
    }
    else {
        pval.points = seq(0.002, 1, .002)
    }
    
    qZ.points = seq(0, 1, .001)
    
    #     weight.matrix = sapply(qZ.points, function(q) dnorm(fq$qZ, q, .01)) / dnorm(0, 0, .01)
    #     pval.matrix = sapply(pval.points, function(p) fq$pvalue < p)
    #     cd = (t(weight.matrix) %*% pval.matrix) / colSums(weight.matrix)
    #     colnames(cd) = pval.points
    #     rownames(cd) = qZ.points
    #     
    #     qval = fq$pi0.func(pval.points) * t(pval.points / t(cd))
    #     qval[qval > 1] = 1
    #     qval.m = melt(qval)
    
    if (!is.null(fq$oracle)) {
        g = g + geom_point(aes(col=oracle), data=fq, size=3)
    }
    else {
        g = g + geom_point(aes(col=confidence),
                           data=fq, size=3)
    }
    # , color=factor(..level..))
    #g = g + stat_contour(aes(x=X2, y=X1, z=value), data=qval.m, breaks=threshold)
    if (is.null(fq$oracle)) {
        g = g + geom_hline(aes(yintercept=pvalue, col=qvalue), data=vlines)
    }
    #scale_color_gradient(name="Q-value", low="red", high="white") +
    
    if (!doublelog) {
        highest.pval = max(fq$pvalue[fq$fqvalue < max(threshold)])
        g = g + ylim(0, highest.pval)
    }

    return(g)
}


#' Compare functional q-value to traditional q-values
#' 
#' Compare traditional qvalue to functional q-values in a scatterplot,
#' colored either by the quantile of Z or by a given oracle
#' 
#' @param fq The data.table returned from fqvalue
#' @param ... additional arguments passed to the traditional qvalue function
#' 
#' @return ggplot object (print to show the plot)
#' 
#' @details If the column \code{oracle} exists in the fqvalue, compareQvalue
#' will use that as the colors of the points. Otherwise, the color depends
#' on the quantiles of Z.
#' 
#' @export
compareQvalue = function(fq, ...) {
    library(ggplot2)
    library(qvalue)
    fq$qvalue = qvalue(fq$pvalue, ...)$qvalue
    g = (ggplot(fq, aes(qvalue, fqvalue,
                        label=round(pvalue, 3))) + xlab("Traditional q-value") +
             ylab("Functional q-value"))
    if (is.null(fq$oracle)) {
        g = g + geom_point(aes(col=qZ))
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
#' 
#' @export
densityCurves = function(fq, cutoff=.05, plottype="lfdr", ...) {
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


#' Plot estimated pi0 as a function of qZ
#' 
#' Plot predicted pi0 as a function of the quantile of Z, and optionally
#' include a horizontal line for the traditional pi0 estimate.
#' 
#' @param fq The data.table returned from fqvalue
#' @param horizontal.line Whether to include a dashed horizontal line
#' representing the traditional estimate of pi0
#' @param ... Extra arguments to give to qvalue for estimating
#' traditional pi0
#' 
#' @return ggplot object (print to show the plot)
#' 
#' @export
plotPi0 = function(fq, horizontal.line=TRUE, ...) {
    library(ggplot2)
    library(qvalue)
    
    pi0 = qvalue(fq$pvalue)$pi0
    g = (ggplot(fq, aes(qZ, pi0)) + geom_line() + xlab("Quantile of Z") +
             ylab(expression(hat(pi[0]))))
    if (horizontal.line) {
        g = g + geom_hline(yintercept=pi0, lty=2)
    }
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


