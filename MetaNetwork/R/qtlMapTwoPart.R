#
# qtlMapTwoPart.R
#
# copyright (c) 2007-2011, GBIC University of Groningen
# last modified mrt, 2010
# first written 2007
# 
# Part of the MetaNetwork Package
#
# Fu J, Swertz MA, Keurentjes JJB, Jansen RC. MetaNetwork: a computational tool for the genetic study of metabolism. Nature Protocols (2007).
#

qtlMapTwoPart <- function (genotypes, traits, spike, filename = NULL){
  sep = ","
  n.marker <- nrow(genotypes)
  if (is.vector(traits)) {
    traits <- t(as.matrix(traits))
  }
  n.traits <- nrow(traits)
  if (!is.null(rownames(traits))) {
    name.traits <- rownames(traits)
  }else{
    name.traits <- paste("trait", 1:n.traits, sep = "")
  }
  if (!is.null(rownames(genotypes))) {
    name.marker <- rownames(genotypes)
  }else{
    name.marker <- paste("M", 1:n.marker, sep = "")
  }
  if (any(colnames(genotypes) != colnames(traits))) {
    stop("Error: mis-matached order of individuals between genotypes and traits")
  }
  if (!is.null(filename)) {
    cat("", rownames(genotypes), file = filename, sep = sep, append = FALSE)
    cat("\n", file = filename, sep = "", append = TRUE)
  }
  lod <- NULL
  for (i in 1:n.traits) {
    qtl.p <- NULL
    for (j in 1:n.marker) {
      z <- traits[i, ]
      z[traits[i, ] <= spike] <- 0
      z[traits[i, ] > spike] <- 1
      za <- z[which(genotypes[j, ] == 1)]
      zb <- z[which(genotypes[j, ] == 2)]
      p1 <- wilcox.test(za, zb)$p.value
      if(is.na(p1)){
        p1 <- 1
      }
      if(p1 == 0){
        p1 <- 2.2e-16
      }
      y <- traits[i, ][traits[i, ] > spike]
      y.genotypes <- genotypes[j, ][traits[i, ] > spike]
      y.temp <- unique(y.genotypes)
      if (length(y.temp[!is.na(y.temp)]) > 1) {
        model <- lm(y ~ y.genotypes)
        p2 <- anova(model)[[5]][1]
      }else{
        p2 <- 1
      }
      if (is.na(p2)){
        p2 <- 1
      }
      if (mean(traits[i, genotypes[j, ] == 1], na.rm = TRUE) > mean(traits[i, genotypes[j, ] == 2], na.rm = TRUE)) {
        qtl.p <- c(qtl.p, log10(p1 * p2))
      }else{
        qtl.p <- c(qtl.p, -1 * log10(p1 * p2))
      }
    }
    if (!is.null(filename)) {
      cat(name.traits[i], qtl.p, file = filename, sep = sep, append = TRUE)
      cat("\n", file = filename, sep = "", append = TRUE)
    }
    lod <- rbind(lod, qtl.p)
  }
  dimnames(lod) <- list(name.traits, name.marker)
  lod
}