#
# qtlCorrZeroOrder.R
#
# copyright (c) 2007-2011, GBIC University of Groningen
# last modified mrt, 2010
# first written 2007
# 
# Part of the MetaNetwork Package
#
# Fu J, Swertz MA, Keurentjes JJB, Jansen RC. MetaNetwork: a computational tool for the genetic study of metabolism. Nature Protocols (2007).
#

qtlCorrZeroOrder <- function (markers, qtlProfiles, qtlThres, filename = NULL){
    QTLcc0 <- function(markers, qtl1, qtl2, qtlThres) {
        CI1 <- qtlSupportInterval(markers, abs(qtl1), qtlThres = qtlThres)
        CI2 <- qtlSupportInterval(markers, abs(qtl2), qtlThres = qtlThres)
        m1 <- NULL
        m2 <- NULL
        if (!is.null(CI1)) {
            for (i in 1:nrow(CI1)) {
                m1 <- c(m1, CI1[i, 1]:CI1[i, 2])
            }
        }
        if (!is.null(CI2)) {
            for (i in 1:nrow(CI2)) {
                m2 <- c(m2, CI2[i, 1]:CI2[i, 2])
            }
            n.CI2 <- nrow(CI2)
        }
        if (is.null(CI1) | is.null(CI2)) {
            r <- 0
        }
        else {
            m <- sort(unique(c(m1, m2)))
            q1 <- q2 <- rep(0, length(qtl1))
            q1[m1] <- qtl1[m1]
            q2[m2] <- qtl2[m2]
            r <- 2 * sum(q1[m] * q2[m])/sum(c(q1[m]^2, q2[m]^2))
        }
        r
    }
    cor.array <- matrix(0, nrow = nrow(qtlProfiles), ncol = nrow(qtlProfiles))
    name <- rownames(qtlProfiles)
    dimnames(cor.array) <- list(name, name)
    cor.array[1, 1] <- 1
    if (!is.null(filename)) {
        cat("", name, file = filename, sep = ",")
        cat("\n", file = filename, append = TRUE)
    }
    for (i in 1:nrow(qtlProfiles)) {
        for (j in 1:i) {
            cor.array[j, i] <- cor.array[i, j] <- QTLcc0(markers = markers, 
                qtl1 = qtlProfiles[i, ], qtl2 = qtlProfiles[j, 
                  ], qtlThres = qtlThres)
        }
    }
    if (!is.null(filename)) {
        for (i in 1:nrow(cor.array)) {
            cat(name[i], cor.array[i, ], file = filename, append = TRUE, sep = ",")
            cat("\n", file = filename, append = TRUE)
        }
    }
    cor.array
}