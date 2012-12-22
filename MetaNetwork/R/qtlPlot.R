#
# qtlPlot.R
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

qtlPlot <- function (markers, qtlProfiles, qtlThres, addTitle = NULL, addMarkerLabels = FALSE, addLegend = TRUE, color = NULL){
    markers.temp <- markers
    n.chr <- length(unique(markers[, 1]))
    if (n.chr > 1) {
        for (i in 2:n.chr) {
            markers.temp[markers[, 1] == i, 2] <- markers.temp[markers[, 
                1] == i, 2] + max(markers.temp[markers[, 1] == 
                i - 1, 2]) + 5
        }
    }
    if (is.matrix(qtlProfiles)) {
        n.trait <- nrow(qtlProfiles)
    }
    else {
        qtlProfiles <- t(as.matrix(qtlProfiles))
        n.trait <- 1
    }
    if (is.null(rownames(qtlProfiles))) {
        rownames(qtlProfiles) <- paste("trait", 1:n.trait, sep = "")
    }
    if (addLegend) {
        x1 <- 0.3 * ceiling(n.trait/3)
    }
    else {
        x1 <- 0
    }
    x2 <- 3.5
    windows(width = 8, height = x2 + x1)
    par(fig = c(0, 1, x1/(x1 + x2), 1))
    par(mai = c(0.8, 0.8, 0.2, 0.2))
    xlim = range(markers.temp[, 2])
    ylim = range(c(qtlProfiles))
    plot(markers.temp[, 2], qtlProfiles[1, ], type = "n", xlim = xlim, ylim = ylim, axes = FALSE, xlab = "Genome", ylab = "-logp QTLs", main = addTitle)
    for (j in 1:1000) {
        lines(c(xlim[1] + (j * xlim[2]/1000), xlim[1] + (j * xlim[2]/1000)), c(ylim[1], ylim[2]), col = grey(0.95), lwd = 5)
    }
    if (is.null(color)) {
        color <- 1:n.trait
    }
    for (i in 1:n.trait) {
        lines(markers.temp[, 2], qtlProfiles[i, ], lwd = 2, col = color[i], 
            lty = i)
    }
    lines(markers.temp[, 2], rep(qtlThres, length = length(markers.temp[, 
        2])), lty = 2)
    lines(markers.temp[, 2], rep(-1 * qtlThres, length = length(markers.temp[, 
        2])), lty = 2)
    lines(markers.temp[, 2], rep(0, length = length(markers.temp[, 
        2])), lty = 1)
    axis(2)
    if (addMarkerLabels) {
        axis(1, markers.temp[, 2], labels = FALSE)
        text(markers.temp[, 2], -0.15 * ylim[2] + ylim[1], labels = rownames(markers), cex = 0.6, adj = 1, srt = 60, xpd = TRUE)
    }else{
        axis(1, markers.temp[, 2], labels = FALSE)
        text(markers.temp[, 2], -0.15 * ylim[2] + ylim[1], labels = round(markers[,2], digits = 1), cex = 0.6, adj = 1, srt = 60, xpd = TRUE)
    }
    if (n.chr > 1) {
        for (k in 1:(n.chr - 1)) {
            x <- mean(c(max(markers.temp[markers.temp[, 1] == k, 2]), min(markers.temp[markers.temp[, 1] == k + 1, 2])))
            lines(c(x, x), c(ylim[1], ylim[2]), col = "white", lwd = 2.5)
            text(x + 20, ylim[2] - (ylim[2] - ylim[1])/20, paste("Chr", k + 1))
        }
    }
    text(20, ylim[2] - (ylim[2] - ylim[1])/20, "Chr 1")
    if (addLegend) {
        par(fig = c(0, 1, 0, x1/(x1 + x2)), new = TRUE)
        par(mai = c(0, 0.8, 0, 0.2))
        plot(c(0, 1), c(0, ceiling(n.trait/3) + 2), type = "n", xlab = "", ylab = "", axes = FALSE)
        if (n.trait == 1) {
            legend(0.5, ceiling(n.trait/3) + 2, legend = rownames(qtlProfiles), col = color, lty = 1:n.trait, bty = "n", xjust = 0.5)
        }
        if (n.trait == 2) {
            legend(0.5, ceiling(n.trait/3) + 2, legend = rownames(qtlProfiles), col = color, lty = 1:n.trait, bty = "n", xjust = 0.5, ncol = 2)
        }
        if (n.trait > 2) {
            y <- 1
            while (y <= n.trait) {
                x <- y%%3
                if (x == 0) {
                  x <- 3
                }
                z <- ceiling(y/3)
                legend((x - 1) * 0.33, ceiling(n.trait/3) - z + 3, legend = rownames(qtlProfiles)[y], col = color[y], lty = y, bty = "n")
                y <- y + 1
            }
        }
    }
}
