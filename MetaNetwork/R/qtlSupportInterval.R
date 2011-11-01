#
# qtlSupportInterval.R
#
# copyright (c) 2007-2011, GBIC University of Groningen
# last modified mrt, 2010
# first written 2007
# 
# Part of the MetaNetwork Package
#
# Fu J, Swertz MA, Keurentjes JJB, Jansen RC. MetaNetwork: a computational tool for the genetic study of metabolism. Nature Protocols (2007).
#

qtlSupportInterval <- function (markers, oneQtlProfile, qtlThres, interval.dropoff = 1.5){
  n.chr <- length(unique(markers[, 1]))
  oneQtlProfile <- abs(oneQtlProfile)
  QTLitv <- NULL
  q.max <- max(oneQtlProfile)
  if (q.max >= qtlThres) {
    for (i in 1:n.chr) {
      chr.mk <- which(markers[, 1] == i)
      chr.qtl <- oneQtlProfile[chr.mk]
      trend <- diff(oneQtlProfile[chr.mk])
      label <- trend
      label[trend < 0] <- -1
      label[trend > 0] <- 1
      turn <- NULL
      for (j in 1:(length(label) - 1)) {
        a <- label[j]
        n <- 0
        while (a == 0 & n < j) {
          n <- n + 1
          a <- label[j - n]
        }
        b <- label[j + 1]
        n <- 1
        while (b == 0 & n < length(label)) {
          n <- n + 1
          b <- label[j + n]
        }
        if (a != 0 & b != 0) {
          turn <- c(turn, a * b)
        }
        if (a == 0 | b == 0) {
          turn <- c(turn, -1)
        }
      }
      turnPoint <- NULL
      if (label[1] <= 0 & chr.qtl[1] >= qtlThres) {
          turnPoint <- c(turnPoint, 1)
      }
      for (k in 1:length(turn)) {
        if (turn[k] < 0 & label[k + 1] < 0 & chr.qtl[k + 1] >= qtlThres) {
          turnPoint <- c(turnPoint, k + 1)
        }
      }
      if (label[k + 1] >= 0 & chr.qtl[k + 2] >= qtlThres) {
        turnPoint <- c(turnPoint, k + 2)
      }
      if (length(turnPoint) > 1) {
        realPeak <- rep(TRUE, length(turnPoint))
        k <- 1
        while (k <= length(which(realPeak))) {
          remain <- turnPoint[realPeak]
          peakValue <- chr.qtl[remain]
          TestWhich <- which(peakValue == sort(peakValue)[k])
          n <- remain[TestWhich]
          if (length(n) > 1) {
            n <- n[1]
          }
          temp.thres <- chr.qtl[n] - interval.dropoff
          t1 <- n
          while (chr.qtl[t1] > temp.thres & t1 > 1) {
            t1 <- t1 - 1
          }
          while (chr.qtl[t1] > temp.thres & t1 > 1) {
            t1 <- t1 - 1
          }
          t2 <- n
          while (chr.qtl[t2] > temp.thres & t2 < length(chr.mk)) {
            t2 <- t2 + 1
          }
          other <- remain[which(1:length(remain) != TestWhich[1])]
          if (any(is.element(other, t1:t2))) {
            realPeak[turnPoint == n] <- FALSE
          }
          else {
            if (chr.qtl[t1] < temp.thres) {
              temp <- chr.mk[t1 + 1]
            }
            else {
              temp <- chr.mk[t1]
            }
            if (chr.qtl[t2] < temp.thres) {
              QTLitv <- rbind(QTLitv, c(temp, chr.mk[t2 - 1]))
            }
            else {
              QTLitv <- rbind(QTLitv, c(temp, chr.mk[t2]))
            }
            k <- k + 1
          }
        }
      }else if (length(turnPoint) == 1) {
        n <- turnPoint
        temp.thres <- chr.qtl[turnPoint] - interval.dropoff
        t1 <- t2 <- n
        while (chr.qtl[t1] > temp.thres & t1 > 1) {
          t1 <- t1 - 1
        }
        while (chr.qtl[t2] > temp.thres & t2 < length(chr.mk)) {
          t2 <- t2 + 1
        }
        if (chr.qtl[t1] < temp.thres) {
          temp <- chr.mk[t1 + 1]
        }
        else {
          temp <- chr.mk[t1]
        }
        if (chr.qtl[t2] < temp.thres) {
          QTLitv <- rbind(QTLitv, c(temp, chr.mk[t2 - 1]))
        }
        else {
          QTLitv <- rbind(QTLitv, c(temp, chr.mk[t2]))
        }
      }
    }
  }
  QTLitv
}