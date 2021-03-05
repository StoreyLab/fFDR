#' One row summary of an analysis by the fFDR methodology
#' 
#' @param object An object of class fqvalue.
#' @param oracle A logical value TRUE/FALSE, where TRUE means
#' the oracle procedure is available, and FALSE the contrary
#' @param ... Extra arguments, not used
#' 
#' @return A tibble with columns summarizing the performance of the fFDR methodology in terms of the fqvalue and FDR, with comparison to those obtained by the traditional qvalue method.
#' 
#' @import qvalue
#' 
#' @export
summary.fqvalue <- function(object, oracle, ...) {
    # work on the fqvalue table
    tab <- object$table
    tab$oracle <- oracle
    
    # if traditional qvalues aren't already present, calculate them
    tab <- tab %>%
        dplyr::mutate(q.value = qvalue(p.value)[["qvalues"]])

    FP <- sum(!tab$oracle[tab$fq.value < .05])
    power <- sum(tab$fq.value < .05)
    test.pval <- binom.test(FP, power, .05, alternative = "greater")$p.value
    
    # summarize within parameters
    tab %>%
        dplyr::summarize(q.value.power = mean(q.value < .05),
                         fq.value.power = mean(fq.value < .05),
                         q.value.fdr = mean(!oracle[q.value < .05]),
                         fq.value.fdr = mean(!oracle[fq.value < .05]),
                         spearman = cor(q.value, fq.value, method = "spearman"),
                         FDR.binom.pval = test.pval)
}
