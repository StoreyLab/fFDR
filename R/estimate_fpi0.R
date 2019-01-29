#' estimate pi0 as a function of z
#' 
#' \code{fpi0} estimates a pi0 function that depends on another
#' parameter z.
#' 
#' @param p a vector of p-values
#' @param z0 values on which pi0 is believed to depend
#' @param lambda Possible choices of lambda, by default {.4, .5, ..., .9}
#' @param method Either "gam" (default), "glm", "kernel", or "bin"
#' @param df Degrees of freedom to use in spline in "gam" method
#' @param breaks Either a number of (evenly spaced) break points in "bin" method,
#' or a vector of break points (from 0 to 1) to use for bins
#' @param ... Additional arguments to glm, gam or the kernel estimator
#' 
#' @details In short, the "glm", "gam", and "kernel" methods attempt
#' to estimate:
#' 
#' \code{pi_0(z) = Pr(p>lambda|z)/(1-lambda)}
#' 
#' The factor z0 that the user provides is transformed to produce z, according to
#' z = rank(z0) / length(z0). This ensures that z has an approximately uniform
#' distribution on the interval (0, 1].
#' 
#' The glm and gam approaches define a variable \code{phi=I{p>lambda}}, and
#' use a modification of logistic regression to fit \code{phi~f(z)}. The
#' kernel density estimate examines the density of z where p>lambda and computes
#' Pr(p>lambda) from that.
#' 
#' Binning simply computes the Storey estimate with the given lambda
#' within each bin.
#' 
#' @return an "fPi0" object, which contains the following components:
#'   \item{table}{a tbl_df with one row for each hypothesis, and columns
#'   \code{p.value}, \code{z}, \code{z0}, and \code{fpi0}}
#'   \item{tableLambda}{An expanded version of the table that shows the
#'   results for each choice of lambda}
#'   \item{MISE}{Calculations on the estimated error of each choice of lambda}
#'   \item{lambda}{Chosen value of lambda}
#' 
#' @importFrom dplyr mutate filter group_by %>%
#' 
#' @export
estimate_fpi0 <- function(p, z0, lambda = seq(.4, .9, .1), method = "gam",
                   df = 3, breaks = 5, ...) {
    # check p-values, and assumptions of model
    if (min(p) < 0 || max(p) > 1) {
        stop("P-values not in valid range")
    }
    if (length(p) != length(z0)) {
        stop("Vectors of p-values and z0 must be same length")
    }

    # transform z0 into quantiles
    z <- rank(z0) / length(z0)
    
    method <- match.arg(method, c("glm", "gam", "kernel", "bin"))

    pi0hat_func <- function(lambda) {
        # set up a function for estimating the fpi0, which
        # returns a list containing fpi0 and the fit object
        phi <- as.numeric(p > lambda)
        fit <- NULL
        
        if (method == "glm") {
            fit <- glm(phi ~ z, family = constrained.binomial(1 - lambda), ...)
            pi0 <- fitted.values(fit) / (1 - lambda)
        } else if (method == "gam") {
            fit <- gam::gam(phi ~ splines::ns(z, df),
                            family = constrained.binomial(1 - lambda), ...)
            pi0 <- fitted.values(fit) / (1 - lambda)
        } else if (method == "kernel") {
            kd <- kernelUnitInterval(z[phi == 1], transformation = "probit",
                                     eval.points = z, ...)
            pi0 <- kd$fx * mean(phi) / (1 - lambda)
        } else if (method == "bin") {
            if (length(breaks) == 1) {
                breaks <- seq(0, 1, 1 / breaks)
            }
            b <- findInterval(z, breaks, rightmost.closed = TRUE)
            freq <- table(b, phi)
            b.pi0 <- (freq[, 2] / rowSums(freq)) / (1 - lambda)
            pi0 <- b.pi0[b]
        }
        pi0
    }
    
    # choose lambda
    fpi0s <- data.frame(p.value = p, z = z) %>%
        tidyr::crossing(lambda = lambda) %>%
        dplyr::group_by(lambda) %>%
        mutate(fpi0 = pi0hat_func(lambda[1]))

    # estimate \hat{phi} using the lowest lambda as reference 
    ref <- fpi0s %>%
        dplyr::ungroup() %>%
        filter(lambda == min(lambda)) %>%
        .$fpi0
    
    # calculate phi.hat, and correct for cases > 1
    fpi0s <- fpi0s %>%
        mutate(k = optimize(function(k) mean((ref - k * (1 - ref) - fpi0) ^ 2),
                            interval = c(-1, 1))$minimum,
               phi.hat = ref - k * (1 - ref),
               fpi0 = pmin(fpi0, 1),
               phi.hat = pmin(phi.hat, 1))

    pi0.S <- qvalue(p)$pi0

    stats <- fpi0s %>%
        dplyr::summarize(omega = mean((fpi0 - phi.hat) ^ 2),
                         delta.sq = (max(mean(fpi0) - pi0.S, 0)) ^ 2) %>%
        mutate(MISE = omega + delta.sq)

    # extract chosen lambda and fpi0 estimate
    lambda.hat = stats$lambda[which.min(stats$MISE)]
    fpi0s <- fpi0s %>%
        mutate(chosen = (lambda == lambda.hat))
    fpi0 <- fpi0s %>%
        filter(chosen == TRUE) %>%
        .$fpi0

    tab = dplyr::data_frame(p.value = p, z = z, z0 = z0, fpi0 = fpi0)
    ret <- list(table = tab, tableLambda = fpi0s, MISE = stats,
                lambda = lambda.hat)
    class(ret) <- "fPi0"
    ret
}


#' Extract functional pi0 estimates.
#' 
#' @param x "fPi0" object
#' @param ... Extra arguments, not used
#' 
#' @export
as.double.fPi0 <- function(x, ...) {
    x$table$fpi0
}


#' Extract functional pi0 estimates.
#' 
#' @param x "fPi0" object
#' @param ... Extra arguments, not used
#' 
#' @export
print.fPi0 <- function(x, ...) {
    cat("Estimation of functional pi0 on", length(x), "p-values",
        "using method", x$method, "with chosen lambda =",
        x$lambda, "\n\n")
    cat("Use plot() on this object to observe how fpi0 varies with z.",
        "Use as.numeric() on this object to access the vector of fpi0",
        "predictions, or as.data.frame() to extract a table comparing",
        "p-values, z, and fpi0.\n")
}



#' constrained binomial family
#' 
#' This is a variation on the binomial family for GLMs and GAMs, where the
#' value of g(x) ranges from 0 to a maximum value (while for the logistic
#' link, it ranges from (0, 1).
#' 
#' @param maximum ceiling constraint on the value of g(x)
#' 
#' @return An object of class "family"
constrained.binomial = function(maximum) {
    link <- structure(list(name = paste0("constrained.logit (0, ", maximum, ")"),
                          linkfun = function(mu) log((mu / maximum) / (1 - mu / maximum)),
                          linkinv = function(eta) maximum / (1 + exp(-eta)),
                          mu.eta = function(eta) maximum / (1 + exp(-eta)),
                          valideta = function(eta) TRUE), class = "link-glm")
    
    fam <- binomial(link)
    fam$validmu <- function(mu) all(mu > 0) && all(mu < maximum)
    fam$family <- paste0("constrained.binomial (0, ", maximum, ")")
    
    # for mgcv
    fam$d2link <- function(mu) { 1 / (1 - (mu / maximum)) ^ 2 - 1 / (mu / maximum) ^ 2 }
    fam$d3link <- function(mu) { 2 / (1 - (mu / maximum)) ^ 3 + 2 / (mu / maximum) ^ 3 }
    fam$d4link <- function(mu) { 6 / (1 - (mu / maximum)) ^ 4 - 6 / (mu / maximum) ^ 4 }
    fam$dvar <- function(mu) rep.int(1, length(mu))
    fam$d3var <- fam$d2var <- function(mu) rep.int(0, length(mu))
    
    # new addition to initialization: mu cannot be greater than maximum
    new.line <- substitute({mustart <- mustart * maximum}, list(maximum = maximum))
    fam$initialize <- c(fam$initialize, new.line)
    fam
}


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
#' @param ... Extra arguments, not used
#' 
#' @return Vector of fPi0 predictions
predict.fPi0 <- function(object, z0 = NULL, z = NULL, lambda = NULL, ...) {
    if (is.null(z) && is.null(z0) && is.null(lambda)) {
        return(object$table$fpi0)
    }
    if (!is.null(z) && !is.null(z0)) {
        stop("Cannot give both z0 and z values to predict FPi0")
    }
    
    if (is.null(lambda)) {
        tab <- object$table
    } else {
        if (!(lambda %in% object$tableLambda$lambda)) {
            stop(paste("Cannot predict with lambda = ", lambda))
        }
        l <- lambda
        tab <- object$tableLambda %>%
            dplyr::filter(lambda == l)
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
}
