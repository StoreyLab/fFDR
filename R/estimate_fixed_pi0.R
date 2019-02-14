#' Estimate a constant functional proportion
#' 
#' Estimate the constant functional proportion pi0(z) when it is independent of the informative variabe z
#' 
#' @param p A vector of p-values
#' @param z0 A vector of observations from the informative variable, of the same length as \code{p}
#' @param lambda Choices of the tuning parameter "lambda" (for the p-value) used to estimate pi0(z)
#' @param tau Choices of the tuning parameter "tau" (for the informative variable) used to estiamte pi0(z)
#' @param pi0.method The method used to estimate pi0(z), being either "smoother" (default) or "bootstrap". The method depends on the tuning parameters "lambda" and "tau"
#' @param lambda.df Degrees of freedom for smoother in argument "lambda"
#' @param tau.df Degrees of freedom for smoother in argument "tau"
#'
#' @details Assume the random variable z0 may affect the power of a statistical test (that induces the p-values) or the likelihood of a true null hypothesis. The m observations z0_i, i=1,...,m of z0 are quantile transformed into z_i, i=1,...,m such that z_i = rank(z0_i) / m, where rank(z0_i) is the rank of z0_i among z0_i, i=1,...,m. Consequently, z_i, i=1,...,m are approximately uniformly distributed on the interval [0,1]. When z_i, i=1,...,m are regarded as observations from the random variable z, then z is approximately uniformly distributed on [0,1]. Namely, z0 has been quantile transformed into z, and they are equivalent. Further, z or z0 is referred to as the informative variable.
#'
#' @return Estimate of pi0(z)
#' 
#' @import splines
#' @importFrom dplyr %>% filter
#'
#' @export
estimate_fixed_pi0 <- function(p, z0, lambda = seq(0, .9, .05), tau = seq(0, .9, .05),
                    pi0.method = "smoother", lambda.df = 3, tau.df = 3) {
    combinations <- fixed_pi0_table(p, z0, lambda, tau, pi0.method,
                                    lambda.df, tau.df)
    
    if (pi0.method == "smoother") {
        pi0 <- combinations %>%
            filter(lambda == max(lambda) & tau == max(tau)) %>%
            .$smoothed
    }
    else if (pi0.method == "bootstrap") {
        pi0 <- combinations$pi0hat[which.min(combinations$MSEhat)]
    }
    return(pi0)
}


#' Provide estimates of the functional proportion at varying values of lambda and tau
#' 
#' Estimate the functional proportion pi0(z) at varying values of lambda and tau when pi0(z) is independent of the informative variable z
#' 
#' @param p A vector of p-values
#' @param z0 A vector of observations from the informative variable, of the same length as \code{p}
#' @param lambda Choices of the tuning parameter "lambda" (for the p-value) used to estimate pi0(z)
#' @param tau Choices of the tuning parameter "tau" (for the informative variable) used to estiamte pi0(z)
#' @param pi0.method The method used to estimate pi0(z), being either "smoother" (default) or "bootstrap". The method depends on the tuning parameters "lambda" and "tau"
#' @param lambda.df Degrees of freedom for smoother in argument "lambda"
#' @param tau.df Degrees of freedom for smoother in argument "tau"
#' 
#' @details Assume the random variable z0 may affect the power of a statistical test (that induces the p-values) or the likelihood of a true null hypothesis. The m observations z0_i, i=1,...,m of z0 are quantile transformed into z_i, i=1,...,m such that z_i = rank(z0_i) / m, where rank(z0_i) is the rank of z0_i among z0_i, i=1,...,m. Consequently, z_i, i=1,...,m are approximately uniformly distributed on the interval [0,1]. When z_i, i=1,...,m are regarded as observations from the random variable z, then z is approximately uniformly distributed on [0,1]. Namely, z0 has been quantile transformed into z, and they are equivalent. Further, z or z0 is referred to as the informative variable.
#' 
#' @importFrom dplyr %>% filter mutate group_by
#' 
#' @return A data.frame with the following columns
#' 
#' \item{lambda}{Choices of lambda as the thresholds for p}
#' \item{tau}{Choices of tau as the thresholds for z}
#' \item{L}{The number of hypotheses for which p > lambda and z > tau}
#' \item{pi0hat}{The estimated pi0(z) using these choices of lambda and tau}
#' \item{variance}{The estimated variance of pi0hat}
#' \item{se}{The estimated standard error of pi0hat}
#' \item{bias}{The estimated bias of pi0hat}
#' \item{MSEhat}{The estimated mean squared error (MSE) of pi0hat}
#' 
#' @export
fixed_pi0_table <- function(p, z0, lambda=seq(0, .9, .05), tau=seq(0, .9, .05),
                         pi0.method=NULL, lambda.df=3, tau.df=3) {
    m <- length(p)
    z <- rank(z0) / m
    combinations <- expand.grid(lambda = lambda, tau = tau)
    
    combinations <- combinations %>%
        group_by(lambda, tau) %>%
        mutate(L = sum(p > lambda & z > tau)) %>%
        group_by() %>%
        mutate(pi0hat = L / (m * (1 - lambda) * (1 - tau)))
    if (is.null(pi0.method) || pi0.method == "smoother") {
        lfit <- lm(pi0hat ~ ns(lambda, lambda.df) + ns(tau, tau.df),
                  combinations)
        combinations$smoothed <- predict(lfit)
    }
    if (is.null(pi0.method) || pi0.method == "bootstrap") {
        combinations <- combinations %>%
            mutate(variance = L / (m * (1 - lambda) * (1 - tau)) ^ 2 * (1 - L / m)) %>%
            mutate(se = sqrt(variance)) %>%
            mutate(bias = pi0hat - quantile(pi0hat, .1)) %>%
            mutate(MSEhat = bias ^ 2 + variance)
    }
    combinations
}
