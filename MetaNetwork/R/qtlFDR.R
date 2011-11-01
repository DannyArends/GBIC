#
# qtlFDR.R
#
# copyright (c) 2007-2011, GBIC University of Groningen
# last modified mrt, 2010
# first written 2007
# 
# Part of the MetaNetwork Package
#
# Fu J, Swertz MA, Keurentjes JJB, Jansen RC. MetaNetwork: a computational tool for the genetic study of metabolism. Nature Protocols (2007).
#

qtlFDR <- function (qtlProfiles, fdrThres = 0.05, qtlThres = NULL){
    pValues <- -1 * as.vector(abs(qtlProfiles))
    pValues <- 10^pValues
    pValues <- as.vector(pValues)
    pThres <- 10^(-1 * qtlThres)
    qValues <- qvalue(pValues)
    result <- NULL
    if (!is.null(fdrThres)) {
        pThres.thres <- max(qValues$pvalue[qValues$qvalue < fdrThres])
        result <- rbind(result, c(fdrThres, pThres.thres, -1 * 
            log10(pThres.thres)))
    }
    if (!is.null(pThres)) {
        maxQValue <- max(qValues$qvalue[qValues$pvalue < pThres])
        result <- rbind(result, c(maxQValue, pThres, -1 * log10(pThres)))
    }
    colnames(result) <- c("qValue(fdr)", "pValue", "-log10P(qtlThres)")
    result
}