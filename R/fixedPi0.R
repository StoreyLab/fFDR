#' compute pi0 where the truth is independent of Z
#' 
#' \code{fixedPi0} estimates pi0 for a set of p-values using a power-informative
#' factor z0. This is the assumption that although power increases with z0,
#' the truth of the hypotheses is independent of it.
#' 
#' @param p a vector of p-values
#' @param z0 a vector of a power-informative factor, of the same length as
#' \code{p}
#' @param lambda Possible values of lambda to try
#' @param tau Possible values of tau to try
#' @param pi0.method Either "smoother" (default) or "bootstrap", the method for
#' choosing pi0 based on the tuning parameters "lambda" and "tau"
#' 
#' @return estimated pi0
#' 
#' @import splines
#' @import dplyr
#'
#' @export
fixedPi0 = function(p, z0, lambda=seq(0, .9, .05), tau=seq(0, .9, .05),
                    pi0.method="smoother", lambda.df=3, tau.df=3) {
    combinations = fixedPi0Table(p, z0, lambda, tau, pi0.method,
                                 lambda.df, tau.df)
    
    if (pi0.method == "smoother") {
        pi0 = (combinations %>% filter(lambda == max(lambda) & tau == max(tau)))$smoothed
    }
    else if (pi0.method == "bootstrap") {
        pi0 = combinations$pi0hat[which.min(combinations$MSEhat)]
    }
    return(pi0)
}


#' compute a table of pi0hat estimates at varying values of lambda and tau
#' 
#' \code{fixedPi0Table} estimates pi0 for a set of p-values using a power-informative
#' factor z0. This is the assumption that although power increases with z0,
#' the truth of the hypotheses is independent of it.
#' 
#' @param p a vector of p-values
#' @param z0 a vector of a power-informative factor, of the same length as
#' \code{p}
#' @param lambda Possible values of lambda to try
#' @param tau Possible values of tau to try
#' 
#' @import data.table
#' @import dplyr
#' 
#' @return a data.frame with the following columns
#' 
#' \item{lambda}{Choices of lambda, the p-value threshold}
#' \item{tau}{Choices of tau, the z threshold}
#' \item{L}{The number of hypotheses for which p > lambda and z > tau}
#' \item{pi0hat}{The estimate of pi0 using these choices of lambda and tau}
#' \item{variance}{The estimated variance of the pi0 estimate}
#' \item{se}{The estimated standard error of the pi0 estimate}
#' \item{bias}{The estimated bias of the pi0 estimate}
#' \item{MSEhat}{The estimated mean squared error of the pi0 estimate, typically used to choose MSE}
#' 
#' @export
fixedPi0Table = function(p, z0, lambda=seq(0, .9, .05), tau=seq(0, .9, .05),
                         pi0.method=NULL, lambda.df=3, tau.df=3) {
    m = length(p)
    z = rank(z0) / m
    combinations = expand.grid(lambda=lambda, tau=tau)
    
    combinations = combinations %>% group_by(lambda, tau) %>%
        mutate(L=sum(p > lambda & z > tau)) %>% group_by() %>%
        mutate(pi0hat=L/(m*(1 - lambda)*(1 - tau)))
    if (is.null(pi0.method) || pi0.method == "smoother") {
        lfit = lm(pi0hat ~ ns(lambda, lambda.df) + ns(tau, tau.df), combinations)
        combinations$smoothed = predict(lfit)
    }
    if (is.null(pi0.method) || pi0.method == "bootstrap") {
        combinations = combinations %>%
            mutate(variance=L/(m*(1 - lambda)*(1 - tau))^2*(1-L/m)) %>%
            mutate(se=sqrt(variance)) %>%
            mutate(bias=pi0hat - quantile(pi0hat, .1)) %>%
            mutate(MSEhat=bias^2+variance)
    }
    combinations
}
