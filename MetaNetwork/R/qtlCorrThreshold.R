#
# qtlCorrThreshold.R
#
# copyright (c) 2007-2011, GBIC University of Groningen
#
# Modified by Danny Arends (2011)
#
# last modified Apr, 2011
# first written 2007
# 
# Part of the MetaNetwork Package
#
# Fu J, Swertz MA, Keurentjes JJB, Jansen RC. MetaNetwork: a computational tool for the genetic study of metabolism. Nature Protocols (2007).
#

qtlCorrThreshold <- function (markers, genotypes, traits, spike, qtlThres, n.permutations = 10000, alpha = 0.05, method = "qtl"){
  n <- ncol(traits)
  traits.name <- colnames(traits)
  corrPermutations <- NULL
  for (i in 1:n.permutations) {
      temp <- NULL
      for (j in 1:nrow(traits)) {
          temp <- rbind(temp, sample(traits[j, ]))
      }
      colnames(temp) <- traits.name
      if (method == "qtl") {
        qtl2 <- qtlMapTwoPart(genotypes = genotypes, traits = temp, spike = spike)
        cor0 <- qtlCorrZeroOrder(markers = markers, qtlProfiles = qtl2, qtlThres = qtlThres)
      }else if (method == "abundance") {
        cor0 <- cor(t(temp), method = "spearman", use = "pairwise.complete.obs")
      }
      cor2 <- qtlCorrSecondOrder(cor0)
      corrPermutations <- c(corrPermutations, max(abs(cor2)))
  }
  corrPermutations <- sort(corrPermutations)
  corrPermutations <<- corrPermutations
  thres <- corrPermutations[length(corrPermutations) * (1 - alpha/(n - 1))]
  thres
}