#' Given an fPi0 object, plot the estimates of pi_0(z) with multiple
#' choices of lambda
#' 
#' @param x fPi0 object
#' @param horizontal.line Whether to draw a horizontal line showing the
#' Storey estimate of pi0 for each lambda (does not depend on z)
#' @param subsample Maximum number of points to plot, otherwise randomly sampled
#' (if too many points are plotted, it is impossible to distinguish dotted
#' and solid lines)
#' @param ... other plotting parameters (not used)
#' 
#' @return a ggplot2 graph plotting pi0 as a function of z
#' 
#' @import ggplot2
#' 
#' @export
plot.fPi0 <- function(x, horizontal.line = FALSE, subsample = 1e4, ...) {
    tab = copy(x$tableLambda)
    tab[, pi0S := mean(pvalue > lambda) / (1 - lambda[1]), by = lambda]
    tab[, Chosen := (x$tableLambda$lambda == x$lambda)]
    
    if (subsample < nrow(tab)) {
        tab <- tab[sample(nrow(tab), subsample)]
    }
    
    g <- ggplot(tab, aes(z, fpi0, col = lambda, group = lambda, lty = Chosen)) +
        geom_line() +
        scale_linetype_manual(values = c(3, 1)) +
        labs(col = expression(lambda)) +
        ylab(expression(hat(pi)[0]^(lambda) ~ (z)))
    
    if (horizontal.line) {
        g <- g + geom_hline(aes(yintercept=pi0S, col=lambda, lty=Chosen))
    }
    return(g)
}


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
#' @param ... Additional arguments (not used)
#' 
#' @import ggplot2
#' @import reshape2
#' @import gridExtra
#' 
#' @export
plot_MISE <- function(x, ...) {
    cm <- melt(x$tableLambda[, list(lambda, z, fpi0, phi.hat)], id=c("lambda", "z"))
    
    # plot comparing fpi0 to the reference at each stage
    cm$lambda <- paste(cm$lambda, ifelse(cm$lambda == x$lambda, "(Chosen)", ""))
    cm$variable <- ifelse(cm$variable == "fpi0", "Estimate", "Reference")
    g1 <- ggplot(cm, aes(z, value, lty = variable)) +
        geom_line() +
        facet_wrap(~ lambda) +
        labs(lty = "") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    # plot comparing the bias, variance and MISE estimates
    sm = melt(x$MISE, id = "lambda")
    g2 = ggplot(sm, aes(lambda, value, col = variable)) +
        geom_line() +
        geom_vline(xintercept = x$lambda, lty = 2) +
        xlab(expression(lambda)) +
        labs(col = "")
    
    return(arrangeGrob(g1, g2, nrow = 2))
}
