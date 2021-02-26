#' Estimate the functional proportion
#' 
#' Estimate the functional proportion pi0(z) when it is not independent of the informative variable z.
#' 
#' @param p A vector of p-values
#' @param z0 A vector of observations from the informative variable, of the same length as \code{p}
#' @param lambda Choices of the tuning parameter "lambda", which are by default {.4, .5, ..., .9}
#' @param method Method for estimating the functional proportion pi0(z), which should be either "gam" (default), "glm", "kernel" or "bin"
#' @param df Degrees of freedom to use for the splines in "gam" method
#' @param breaks Either a number of (evenly spaced) break points for "bin" method,
#' or a vector of break points (from 0 to 1) to use for bins
#' @param maxit Max iterations to perform when using the "glm" method.
#' @param ... Additional arguments to "glm", "gam" or "kernel" method
#' 
#' @details Assume the random variable z0 may affect the power of a statistical test (that induces the p-values) or the likelihood of a true null hypothesis. The m observations z0_i, i=1,...,m of z0 are quantile transformed into z_i, i=1,...,m such that z_i = rank(z0_i) / m, where rank(z0_i) is the rank of z0_i among z0_i, i=1,...,m. Consequently, z_i, i=1,...,m are approximately uniformly distributed on the interval [0,1]. When z_i, i=1,...,m are regarded as observations from the random variable z, then z is approximately uniformly distributed on [0,1]. Namely, z0 has been quantile transformed into z, and they are equivalent. Further, z or z0 is referred to as the informative variable.
#' 
#' In short, the "glm", "gam", and "kernel" methods attempt
#' to estimate:
#' 
#' \code{pi_0(z) = Pr(p>lambda|z)/(1-lambda)}
#' 
#' The "glm" and "gam" approaches define an indicate variable \code{phi=I{p>lambda}}, and
#' use a modification of logistic regression to fit \code{phi~f(z)}. The
#' "kernel" method examines the density of z where p>lambda and computes
#' Pr(p>lambda) from that.
#' 
#' Binning simply computes the Storey pi0 estimate with the given lambda
#' within each bin.
#' 
#' @return an "fPi0" object, which contains the following components:
#'   \item{table}{a tibble with one row for each hypothesis, and columns
#'   \code{p.value}, \code{z}, \code{z0}, and \code{fpi0}}
#'   \item{tableLambda}{An expanded version of the table that shows the estimation
#'   results for each choice of lambda}
#'   \item{MISE}{The estimated mean integrated squared error (MISE) of the estimated pi0(z) for each choice of lambda}
#'   \item{lambda}{The chosen value of lambda}
#' 
#' @importFrom dplyr mutate filter group_by %>%
#' 
#' @examples 
#' sim.ttests = simulate_t_tests(m = 1000)
#' fpi0 <- estimate_fpi0(p = sim.ttests$p.value, z0 = sim.ttests$n, method = "kernel")
#' 
#' @export
estimate_fpi0 <- function(p, z0, lambda = seq(.4, .9, .1), method = "gam",
                   df = 3, breaks = 5, maxit = 1000, ...) {
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
            fit <- glm(phi ~ z, family = constrained.binomial(1 - lambda), maxit = maxit, ...)
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
        # add temp row index in case (p,z) tuples are non-unique
        dplyr::mutate(i = 1:n()) %>%
        tidyr::crossing(lambda = lambda) %>%
        dplyr::select(-i) %>%
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

    tab = tibble::tibble(p.value = p, z = z, z0 = z0, fpi0 = fpi0)
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
    cat("Estimation of functional proportion on", length(x), "p-values",
        "using method", x$method, "with chosen lambda =",
        x$lambda, "\n\n")
    cat("Use plot() on this object to observe how the estimated pi0(z) varies with z.",
        "Use as.numeric() on this object to access the vector of estimated pi0(z)",
        "predictions, or as.data.frame() to extract a table comparing",
        "p-values, z, and estimated pi0(z).\n")
}



#' constrained binomial family
#' 
#' This is a variation on the binomial family for GLMs and GAMs, where the
#' value of g(x) ranges from 0 to a maximum value (while for the logistic
#' link, it ranges from (0, 1).
#' 
#' @param maximum Ceiling constraint on the value of g(x)
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
    new.line <- substitute(mustart <- mustart * maximum, list(maximum = maximum))
    # convert call to expressiion so it be initialized by eval
    fam$initialize <- as.expression(c(fam$initialize, new.line))

    fam
}


#' Predict the functional proportion
#' 
#' Predict the values of the functional proportion when it is evaluated at a vector z' whose entries are between 0 and 1.
#' 
#' @param object fPi0 object
#' @param z0 New z0 values (when they are between 0 and 1)
#' @param z New z values (as the quantile transformed values of z0)
#' @param lambda The value of the tuning parameter to be used for prediction. If null, defaults to the lambda
#' chosen in the fPi0 object
#' @param ... Extra arguments, not used
#' 
#' @details Assume the random variable z0 may affect the power of a statistical test (that induces the p-values) or the likelihood of a true null hypothesis. The m observations z0_i, i=1,...,m of z0 are quantile transformed into z_i, i=1,...,m such that z_i = rank(z0_i) / m, where rank(z0_i) is the rank of z0_i among z0_i, i=1,...,m. Consequently, z_i, i=1,...,m are approximately uniformly distributed on the interval [0,1]. When z_i, i=1,...,m are regarded as observations from the random variable z, then z is approximately uniformly distributed on [0,1]. Namely, z0 has been quantile transformed into z, and they are equivalent. Further, z or z0 is referred to as the informative variable.
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
