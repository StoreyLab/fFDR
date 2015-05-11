#' One row summary of an fqvalue analysis
#' 
#' @param object An object of class fqvalue.
#' @param oracle A logical vector giving the oracle, where TRUE means
#' the alternative hypothesis holds and FALSE means the null does
#' @param ... Extra arguments, not used
#' 
#' @return A data.table with columns summarizing the performance of fqvalue
#' 
#' @import data.table
#' @import qvalue
#' 
#' @export
summary.fqvalue <- function(object, oracle, ...) {
    # work on the fqvalue table
    tab = object$table
    tab$oracle = oracle
    
    # if traditional qvalues aren't already present, calculate them
    tab[, q.value := qvalue(p.value)$qvalues]
    
    FP = sum(!tab$oracle[tab$fq.value < .05])
    power = sum(tab$fq.value < .05)
    test.pval = binom.test(FP, power, .05, alternative = "greater")$p.value
    
    # summarize within parameters
    tab[, list(q.value.power = mean(q.value < .05),
               fq.value.power = mean(fq.value < .05),
               q.value.fdr = mean(!oracle[q.value < .05]),
               fq.value.fdr = mean(!oracle[fq.value < .05]),
               spearman = cor(q.value, fq.value, method="spearman"),
               FDR.binom.pval = test.pval)]
}
