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
