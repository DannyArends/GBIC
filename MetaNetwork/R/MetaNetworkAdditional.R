#
# MetaNetworkAdditional.R
#
# copyright (c) 2007-2011, GBIC University of Groningen
#
# Modified by Danny Arends (2011)
#
# last modified Apr, 2011
# first written 2011
# 
# Part of the MetaNetwork Package
# Contains: loadData, timeToOutfile, qtlThreshold
#
# Fu J, Swertz MA, Keurentjes JJB, Jansen RC. MetaNetwork: a computational tool for the genetic study of metabolism. Nature Protocols (2007).
#

loadData <- function(filename){
  obj <- read.csv(filename)
  rownames(obj) <- obj[, 1]
  obj <- as.matrix(obj[, -1])
  obj
}


timeToOutfile <- function(t1, t2, outputdir, outfile){
  t2 <- proc.time()
  if ((t2 - t1)[3] > 60){
    cat("\t", "process time", round((t2 - t1)[3]/60, digits = 2), "min\n\n")
    if (!is.null(outputdir)) cat("\t", "process time", round((t2 - t1)[3]/60, digits = 2), "min\n\n", file = outfile, append = TRUE)
  }else{
    cat("\t", "process time", (t2 - t1)[3], "sec\n\n")
    if (!is.null(outputdir)) cat("\t", "process time", (t2 - t1)[3], "sec\n\n", file = outfile, append = TRUE)
  }
  flush.console()
}

qtlThreshold <- function (genotypes, traits, spike, n.simulations = 1000, alpha = 0.05){
  spikeP <- function(x) { length(which(x <= spike)) }
  sdP2 <- function(x) { sd(x[x > spike], na.rm = TRUE) }
  mP2 <- function(x) { mean(x[x > spike], na.rm = TRUE) }
  
  n1 <- apply(traits, 1, spikeP)
  sd2 <- apply(traits, 1, sdP2)
  m2 <- apply(traits, 1, mP2)
  simulations <- NULL
  for (i in 1:n.simulations) {
    values <- sample(c(rnorm(ncol(traits) - ceiling(median(n1)), median(m2), median(sd2)), rep(spike, ceiling(median(n1)))))
    result <- qtlMapTwoPart(genotypes = genotypes, traits = values, spike = spike)
    simulations <- c(simulations, max(abs(result)))
  }
  thres <- sort(simulations)[(1 - alpha) * n.simulations]
  thres
}