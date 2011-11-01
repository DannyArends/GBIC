#
# qtlSummary.R
#
# copyright (c) 2007-2011, GBIC University of Groningen
# last modified mrt, 2010
# first written 2007
# 
# Part of the MetaNetwork Package
#
# Fu J, Swertz MA, Keurentjes JJB, Jansen RC. MetaNetwork: a computational tool for the genetic study of metabolism. Nature Protocols (2007).
#

qtlSummary <- function (markers, genotypes, traits, qtlProfiles, spike, qtlThres, interval.dropoff = 1.5, filename = NULL){
  n.marker <- nrow(genotypes)
  if (is.vector(traits)) {
    traits <- t(as.matrix(traits))
  }
  n.traits <- nrow(traits)
  name.traits <- rownames(traits)
  if (is.null(name.traits)) {
    name.traits <- 1:n.traits
  }
  if (!is.null(filename)) {
    cat("traitName", paste("QTLchr", "(", qtlThres, ")", sep = ""), "QTLmk", "QTLleft(cm)", "QTLpeak(cm)", paste("QTLright(cm)", "(", interval.dropoff, ")", sep = ""), "logp", "VarP1(%)", "VarP2(%)", "Additive", "\n", file = filename, sep = ",")
  }
  out <- NULL
  for (i in 1:n.traits) {
    qtl <- abs(qtlProfiles[i, ])
    itv <- qtlSupportInterval(markers, qtl, qtlThres = qtlThres, interval.dropoff = interval.dropoff)
    z <- traits[i, ]
    z[traits[i, ] > spike] <- 1
    z[traits[i, ] <= spike] <- 0
    if (!is.null(itv)){
      for (j in 1:nrow(itv)) {
        var1 <- 0
        var2 <- 0
        mk1 <- itv[j, 1]
        mk2 <- itv[j, 2]
        chr <- markers[mk1, 1]
        cm1 <- markers[mk1, 2]
        cm2 <- markers[mk2, 2]
        chr.mk <- range(which(markers[, 1] == chr))
        if (mk1 > chr.mk[1]) {
          cm1 = mean(c(markers[mk1, 2], markers[mk1 - 1, 2]))
        }
        if (mk2 < chr.mk[2]) {
          cm2 <- mean(c(markers[mk2, 2], markers[mk2 + 1, 2]))
        }
        logP <- max(qtl[mk1:mk2])
        n <- which(qtl[mk1:mk2] == logP) + mk1 - 1
        n <- n[1]
        mk <- n
        mk.cm <- markers[n, 2]
        if (length(which(z == 0)) > 0) {
          temp <- genotypes[n, ]
          model <- glm(z ~ temp, na.action = "na.omit", family = binomial(link = "probit"))
          var1 <- anova(model)[[2]][2]/anova(model)[[4]][1]
        }
        y <- traits[i, ][traits[i, ] > spike]
        y.genotypes <- genotypes[n, ][traits[i, ] > spike]
        y.temp <- unique(y.genotypes)
        if (length(y.temp[!is.na(y.temp)]) > 1) {
          model2 <- lm(y ~ y.genotypes)
          result <- anova(model2)
          var2 <- result[[2]][1]/(result[[2]][1] + result[[2]][2])
        }
        additive <- 0.5 * (mean(traits[i, genotypes[n, ] == 2], na.rm = TRUE) - mean(traits[i, genotypes[n, ] == 1], na.rm = TRUE))
        out <- rbind(out, data.frame(traitName = name.traits[i], QTLchr = chr, QTLmk = rownames(markers)[mk], QTLleftcm = cm1, QTLpeakcm = markers[mk, 2], QTLrightcm = cm2, logP = round(logP, digits = 1), VarP1 = round(100 * var1, digits = 1), VarP2 = round(100 * var2, digits = 1), additive = round(additive, digits = 1)))
        if (!is.null(filename)) {
          cat(name.traits[i], chr, rownames(markers)[mk], cm1, mk.cm, cm2, round(logP, digits = 1), round(100 * var1, digits = 1), round(100 * var2, digits = 1), round(additive, digits = 1), "\n", file = filename, sep = ",", append = TRUE)
        }
      }
    }else{
      if (!is.null(filename)) {
        cat(name.traits[i], "NS", "\n", file = filename, sep = ",", append = TRUE)
      }
    }
  }
  out
}