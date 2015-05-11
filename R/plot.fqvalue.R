#' log(-log(x)) transformation for a ggplot2 scale
#' 
#' @description this code is a modified version of this:
#' http://stackoverflow.com/questions/11053899
#' 
#' @param base Base used in the exponentiation
reverseloglog_trans <- function(base = exp(1)) {
    trans <- function(x) log(-log(x, base), base)
    inv <- function(x) base ^ (-(base ^ x))
    trans_new(paste0("reverseloglog-", format(base)), trans, inv, 
              function(x) c(.9985, .99, .93, .65, .1, 1e-05,
                            1e-25, 1e-100), 
              domain = c(1e-300, Inf))
}


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
#' @param ... Extra arguments (not used)
#' 
#' @return ggplot object (print to show the plot)
#' 
#' @import ggplot2
#' @import scales
#' @import reshape2
#' @import gridExtra
#' 
#' @export
plot.fqvalue <- function(x, pi0 = TRUE, threshold = c(.005, .01, .05, .1),
                   cloglog = FALSE, ...) {
    tab <- x$table
    qv <- qvalue(tab$p.value)$qvalue
    
    num.below <- sapply(threshold, function(q) sum(qv < q))
    vlines <- data.frame(q.value = as.factor(threshold),
                         p.value = c(0, sort(tab$p.value))[num.below + 1])
    qp <- factor(c(threshold, 1))
    tab$confidence <- qp[findInterval(tab$fq.value, threshold) + 1]
    
    if (!any(tab$fq.value < max(threshold))) {
        # if nothing is significant, just show everything
        highest.pval <- 1
    } else {
        highest.pval <- max(tab$p.value[tab$fq.value < max(threshold)])
    }
    tab <- tab[tab$p.value <= max(highest.pval, max(vlines$p.value)), ]
    
    g <- ggplot(tab, aes(x = z, y = p.value)) +
        scale_color_brewer(palette = "Spectral") +
        theme(axis.text.x = element_text(angle = 30, hjust = 1))
    
    if (cloglog) {
        g <- g + scale_y_continuous(trans = reverseloglog_trans(10))
        pval.points <- 10 ^ (-10 ^ (seq(-4.5, 2, .05)))
    }
    else {
        pval.points <- seq(0.002, 1, .002)
    }
    
    qZ.points <- seq(0, 1, .001)
    
    colorer <- ifelse(is.null(tab$oracle), "confidence", "oracle")
    g <- g + geom_point(aes_string(col = colorer), data = tab, size = 3)

    if (is.null(tab$oracle)) {
        g <- g + geom_hline(aes(yintercept = p.value, col = q.value), data = vlines)
    }

    if (pi0) {
        return(arrangeGrob(g, plot(x$fPi0), nrow = 2))
    }
    
    return(g)
}


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
compare_qvalue = function(fq, ...) {
    tab <- fq$table
    tab$q.value <- qvalue(tab$p.value, ...)$qvalue
    g <- ggplot(tab, aes(q.value, fq.value, label = round(p.value, 3))) +
        geom_abline(col = "red") +
        xlab("Traditional q-value") +
        ylab("Functional q-value")
    
    colorer <- ifelse(is.null(tab$oracle), "z", "oracle")
    g <- g + geom_point(aes_string(col = colorer))
    g
}
