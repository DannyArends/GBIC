#
# findPeakMultiplicity.R
#
# copyright (c) 2007-2011, GBIC University of Groningen
#
# Modified by Danny Arends (2011)
#
# last modified Apr, 2011
# first written 2011
# 
# Part of the MetaNetwork Package
#
# Fu J, Swertz MA, Keurentjes JJB, Jansen RC. MetaNetwork: a computational tool for the genetic study of metabolism. Nature Protocols (2007).
#

findPeakMultiplicity <- function(corrZeroOrder, peaks, corrThres = 0.95, filename = NULL){
  grouping <- function(pairPeaks, starts) {
    anyLink <- NULL
    for (i in 1:nrow(pairPeaks)) {
      anyLink <- c(anyLink, any(is.element(pairPeaks[i, ], starts)))
    }
    anyLink
  }
  isotope <- c(1, 2)
  diffCharged <- c(0.33, 0.5, 2, 3)
  peakMultiplicity <- NULL
  name <- rownames(corrZeroOrder)
  for (i in 2:nrow(corrZeroOrder)) {
      for (j in 1:(i - 1)) {
          m1 <- peaks[which(rownames(peaks) == name[i]), 1]
          m2 <- peaks[which(rownames(peaks) == name[j]), 1]
          massRatio <- round(m1/m2, digits = 2)
          massDiff <- m1 - m2
          relationship <- ""
          if (corrZeroOrder[i, j] > corrThres) {
              if (is.element(massRatio, diffCharged)) {
                relationship <- "diffCharged"
              }
              if (is.element(abs(massDiff), isotope)) {
                relationship <- "isotope"
              }
              if (relationship != "") {
                peakMultiplicity <- rbind(peakMultiplicity, data.frame(peak1 = name[i], mass_charge1 = m1, peak2 = name[j], mass_charge2 = m2, corrCoef = corrZeroOrder[i,j], massDiff = massDiff, massRatio = massRatio, relationship = relationship))
              }
          }
      }
  }
  pairPeaks <- as.matrix(peakMultiplicity[, 1:2])
  clusterNo <- rep(0, nrow(pairPeaks))
  j <- 0
  while (any(clusterNo == 0)) {
      j <- j + 1
      n <- which(clusterNo == 0)[1]
      tempMember <- pairPeaks[n, 1:2]
      anyLink <- grouping(pairPeaks, tempMember)
      m <- which(anyLink)
      while (length(m) > length(n)) {
          tempMember <- unique(as.vector(pairPeaks[m, ]))
          n <- m
          anyLink <- grouping(pairPeaks, tempMember)
          m <- which(anyLink)
      }
      clusterNo[n] <- j
  }
  peakMultiplicity <- cbind(cluster = clusterNo, peakMultiplicity)
  peakMultiplicity <- peakMultiplicity[order(peakMultiplicity[,1]), ]
  if (!is.null(filename)) {
      write.table(peakMultiplicity, file = filename, sep = ",", row.names = FALSE)
  }
  peakMultiplicity
}
