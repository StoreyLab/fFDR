#' Plot the estimated functional proportion
#'
#' Plot the estimated functional proportion for different choices of tuning parameter "lambda".
#'  
#' @param x fPi0 object
#' @param horizontal.line A logical variable TRUE/FALSE, indicating whether to draw a horizontal line showing the
#' Storey estimate of a constant pi0(z) for each choice of lambda
#' @param subsample Maximum number of points to plot; otherwise randomly sampled
#' (if too many points are plotted, it is impossible to distinguish dotted
#' and solid lines)
#' @param ... Other plotting parameters (not used)
#' 
#' @return A ggplot2 graph plotting the functional proportion pi0(z) as a function of the informative variable z. For information on z, please refer to \code{fqvalue} or \code{estimate_fpi0}
#' 
#' @import ggplot2
#' @importFrom dplyr %>%
#' 
#' @export
plot.fPi0 <- function(x, horizontal.line = FALSE, subsample = 1e4, ...) {
    tab <- x$tableLambda %>%
        dplyr::group_by(lambda) %>%
        dplyr::mutate(pi0S = mean(p.value > lambda) / (1 - lambda[1])) %>%
        dplyr::mutate(Chosen = (lambda == x$lambda)) %>%
        dplyr::ungroup()

    if (subsample < nrow(tab)) {
        tab <- tab %>% dplyr::sample_n(subsample)
    }
    
    g <- ggplot(tab, aes(z, fpi0, col = lambda, group = lambda, lty = Chosen)) +
        geom_line() +
        scale_linetype_manual(values = c(3, 1)) +
        labs(col = expression(lambda)) +
        ylab(expression(paste(hat(pi)[0],"(",z, ";",lambda,")",sep="")))
    
    if (horizontal.line) {
        g <- g + geom_hline(aes(yintercept = pi0S, col = lambda, lty = Chosen))
    }
    return(g)
}


#' Plot demonstrating how the tuning parameter "lambda" is determined
#' 
#' Construct a plot containing two panels. The first compares the
#' estimate of the functional proportion pi0(z) for each choice of lambda to a "Reference" estimate of pi0(z) obtained with
#' the smallest of the choices of lambda value. The Reference estimate usually accurately captures the shapce of pi0(z).
#' The second panel shows how omega (the estimated integrated variance of the estimated pi0(z)), delta.sq
#' (the square of the estimated integrated bias of the estimated pi0(z)), and their sum (the estimated mean integrated squared error (MISE)) change with changing
#' lambda. The software's choice of lambda is highlighted as a vertical
#' dashed line.
#' 
#' @param x fpi0 object
#' @param ... Additional arguments (not used)
#' 
#' @import ggplot2
#' @import gridExtra
#' @importFrom dplyr %>%
#' 
#' @export
plot_MISE <- function(x, ...) {
    cm <- x$tableLambda %>%
        dplyr::select(lambda, z, fpi0, phi.hat) %>%
        reshape2::melt(id = c("lambda", "z"))

    # plot comparing fpi0 to the reference at each stage
    cm$lambda <- paste(cm$lambda, ifelse(cm$lambda == x$lambda, "(Chosen)", ""))
    cm$variable <- ifelse(cm$variable == "fpi0", "Estimate", "Reference")
    g1 <- ggplot(cm, aes(z, value, lty = variable)) +
        geom_line() +
        facet_wrap(~ lambda) +
        labs(lty = "") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    # plot comparing the bias, variance and MISE estimates
    sm <- melt(x$MISE, id = "lambda")
    g2 <- ggplot(sm, aes(lambda, value, col = variable)) +
        geom_line() +
        geom_vline(xintercept = x$lambda, lty = 2) +
        xlab(expression(lambda)) +
        labs(col = "")

    # grid::grid.draw(arrangeGrob(g1, g2, nrow = 2))
    return(arrangeGrob(g1, g2, nrow = 2))
    
}
