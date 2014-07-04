# methods for an fPi0 object (besides plots, which are handled in plots.R)

setMethod("length", "fPi0", function(x) nrow(x@table))
setMethod("as.numeric", "fPi0", function(x, ...) x@table$fpi0)
setMethod("as.data.frame", "fPi0", function(x, ...) as.data.frame(x@table))

setMethod("show", "fPi0",
    function(object) {
        cat("Estimation of functional pi0 on", length(object), "pvalues",
            "using method", object@method, "with chosen lambda =",
            object@lambda, "\n\n")
        cat("Use plot() on this object to observe how fpi0 varies with z.",
            "Use as.numeric() on this object to access the vector of fpi0",
            "predictions, or as.data.frame() to extract a table comparing",
            "p-values, z, and fpi0.\n")
    })

#' Predict pi0 for a given z value
#' 
#' Can be given either z0 values (same scale as was given to estFPi0)
#' or z values (after z0 was transformed to uniform; this is the scale
#' on which plots of fpi0 are made)
#' 
#' @param object fPi0 object
#' @param z0 new z0 values (before transforming to uniform)
#' @param z new z values (after transforming to uniform)
#' @param lambda Lambda used for prediction. If null, defaults to the lambda
#' chosen in the fPi0 object
#' 
#' @return Vector of fPi0 predictions
setMethod("predict", "fPi0",
    function(object, z0=NULL, z=NULL, lambda=NULL, ...) {
        if (is.null(z) && is.null(z0) && is.null(lambda)) {
            return(object@table$fpi0)
        }
        if (!is.null(z) && !is.null(z0)) {
            stop("Cannot give both z0 and z values to predict FPi0")
        }

        if (is.null(lambda)) {
            tab = object@table
        } else {
            if (!(lambda %in% object@tableLambda$lambda)) {
                stop(paste("Cannot predict with lambda = ", lambda))
            }
            l = lambda
            tab = object@tableLambda[lambda == l, ]
        }

        # approximate with linear interpolation based on the table
        if (!is.null(z0)) {
            return(approx(tab$z0, tab$fpi0, z0))
        } else if (!is.null(z)) {
            return(approx(tab$z, tab$fpi0, z))
        }
        else {
            return(tab$fpi0)
        }
    })

## Plots

#' Given an fPi0 object, plot the estimates of pi_0(z) with multiple
#' choices of lambda
#' 
#' @param x fPi0 object
#' @param horizontal.line Whether to draw a horizontal line showing the
#' Storey estimate of pi0 for each lambda (does not depend on z)
#' @param ... other plotting parameters (not used)
#' 
#' @return a ggplot2 graph plotting pi0 as a function of z
#' 
#' @import ggplot2
#' 
#' @export
setMethod("plot", c("fPi0"),
          function(x, horizontal.line=FALSE, ...) {
              tab = copy(x@tableLambda)
              tab[, pi0S:=mean(pvalue > lambda) / (1 - lambda[1]), by=lambda]
              tab[, Chosen:=(x@tableLambda$lambda == x@lambda)]
              
              g = ggplot(tab, aes(z, fpi0, col=lambda, group=lambda,
                                  lty=Chosen)) + geom_line()
              g = g + scale_linetype_manual(values=c(3, 1))
              g = g + labs(col=expression(lambda))
              g = g + ylab(expression(hat(pi)[0]^(lambda) ~ (z)))
              
              if (horizontal.line) {
                  g = g + geom_hline(aes(yintercept=pi0S, col=lambda, lty=Chosen))
              }
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
    cm = melt(x@tableLambda[, list(lambda, z, fpi0, phi.hat)], id=c("lambda", "z"))
    
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
    
    return(arrangeGrob(g1, g2, nrow=2))
})
