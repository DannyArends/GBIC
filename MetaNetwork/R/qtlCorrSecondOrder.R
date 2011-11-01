#
# qtlCorrSecondOrder.R
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

qtlCorrSecondOrder <- function (corrZeroOrder, topCorNo = 20, filename = NULL) 
{
    if (!is.null(filename)) {
        cat("", rownames(corrZeroOrder), file = filename, append = FALSE, sep = ",")
        cat("\n", file = filename, append = TRUE)
    }
    secondCor <- matrix(0, nrow = nrow(corrZeroOrder), ncol = ncol(corrZeroOrder))
    dimNo <- nrow(corrZeroOrder)
    n <- 1:dimNo
    for (i in 2:nrow(corrZeroOrder)) {
        for (j in 1:(i - 1)) {
            rxy <- 1
            if (dimNo > topCorNo) {
                avgCor <- apply(abs(corrZeroOrder[c(i, j), ]), 
                  2, mean)
                n <- order(avgCor, decreasing = TRUE)[1:topCorNo]
            }
            m <- 0
            for (x in 2:length(n)) {
                k <- n[x]
                for (y in 1:(x - 1)) {
                  p <- n[y]
                  m <- m + 1
                  if (!is.element(k, c(i, j)) & !is.element(p, 
                    c(i, j)) & corrZeroOrder[i, k] != 1 & corrZeroOrder[j, 
                    k] != 1 & corrZeroOrder[p, k] != 1) {
                    rxyk <- (corrZeroOrder[i, j] - corrZeroOrder[i, 
                      k] * corrZeroOrder[j, k])/sqrt((1 - corrZeroOrder[i, 
                      k]^2) * (1 - corrZeroOrder[j, k]^2))
                    rxpk <- (corrZeroOrder[i, p] - corrZeroOrder[i, 
                      k] * corrZeroOrder[p, k])/sqrt((1 - corrZeroOrder[i, 
                      k]^2) * (1 - corrZeroOrder[p, k]^2))
                    rypk <- (corrZeroOrder[j, p] - corrZeroOrder[j, 
                      k] * corrZeroOrder[p, k])/sqrt((1 - corrZeroOrder[j, 
                      k]^2) * (1 - corrZeroOrder[p, k]^2))
                    if (abs(rxyk) < 1 & abs(rxpk) < 1 & abs(rypk) < 
                      1) {
                      rxykp <- (rxyk - rxpk * rypk)/sqrt((1 - 
                        rxpk^2) * (1 - rypk^2))
                      if (abs(rxy) > abs(rxykp)) {
                        rxy <- rxykp
                      }
                    }
                  }
                }
            }
            secondCor[i, j] <- secondCor[j, i] <- rxy
        }
    }
    dimnames(secondCor) <- dimnames(corrZeroOrder)
    if (!is.null(filename)) {
        for (i in 1:nrow(secondCor)) {
            cat(rownames(secondCor)[i], secondCor[i, ], file = filename, append = TRUE, sep = ",")
            cat("\n", file = filename, append = TRUE)
        }
    }
    secondCor
}