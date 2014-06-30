#' compute pi0 where the truth is independent of Z
#' 
#' \code{fixedPi0} estimates pi0 for a set of p-values using a power-informative
#' factor Z. This  the assumption that although power depends on Z,
#' the truth of the hypotheses is independent of it.
#' 
#' @param p a vector of p-values
#' @param Z a vector of a power-informative factor, of the same length as
#' \code{p}
#' @param lambda Possible values of lambda to try
#' @param tau Possible values of tau to try
#' 
#' @import data.table
#' 
#' @return estimated pi0
#' 
#' @export
fixedPi0 = function(p, Z, lambda=seq(0, .9, .05),
                    tau=c(0, .9, .05)) {
    m = length(p)
    qZ = rank(Z) / m
    combinations = as.data.table(expand.grid(lambda=lambda, tau=tau))
    combinations[, L:=sum(p > lambda & qZ > tau), by=c("lambda", "tau")]
    combinations[, pi0hat:=L/(m*(1 - lambda)*(1 - tau))]
    combinations[, variance:=L/(m*(1 - lambda)*(1 - tau))^2*(1-L/m)]
    combinations[, se:=sqrt(variance)]
    combinations[, pi0hat.1se:=pi0hat+se]

    combinations[, bias:=pi0hat - min(pi0hat)]
    #combinations[, bias:=ifelse(pi0hat < min(pi0hat.1se), 0, pi0hat - min(pi0hat.1se))]
    combinations[, MSEhat:=bias^2+variance]
}
