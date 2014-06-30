#' Summarize results of an fqvalue analysis
#' 
#' @param fq An object of class fqvalue.
#' 
#' @return A data.table with columns summarizing the performance of fqvalue
#' 
#' @import data.table
#' @import qvalue
#' 
#' @export
summary.fqvalue = function(fq) {
    # there may be parameters, allowing for multiple fqvalue objects. These would be the first
    # columns of fq, before pvalue.
    pval.i = which(colnames(fq) == "pvalue")[1]
    if (pval.i == 1) {
        pars = NULL
    } else {
        pars = colnames(fq)[1:(pval.i - 1)]
    }

    # if traditional qvalues aren't already present, calculate them
    fq[, qvalue:=qvalue(pvalue)$qvalue, by=pars]

    FP = sum(!fq$oracle[fq$fqvalue < .05])
    power = sum(fq$fqvalue < .05)
    test.pval = binom.test(FP, power, .05, alternative="greater")$p.value

    # summarize within parameters
    fq[, list(qvalue.power=mean(qvalue < .05), fqvalue.power=mean(fqvalue < .05),
              qvalue.fdr=mean(!oracle[qvalue < .05]),
              fqvalue.fdr=mean(!oracle[fqvalue < .05]),
              spearman=cor(qvalue, fqvalue, method="spearman"),
              FDR.binom.pval=test.pval), by=pars]
}
