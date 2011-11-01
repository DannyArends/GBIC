#
# MetaNetwork.R
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

MetaNetwork <- function(markers, genotypes, traits, spike, qtlProfiles = NULL, 
                        peaks = NULL, qtlThres = NULL, qtlSumm = NULL, corrZeroOrder = NULL, 
                        corrSecondOrder = NULL, corrMethod = "qtl", corrThres = 0, 
                        cytoFiles = TRUE, outputdir = "./MetaNetwork"){
  Mdif <- apply(markers, 2, diff)
  if (any(Mdif < 0)) {
    chrDif <- diff(markers[, 1])
    if (any(!is.element(chrDif, c(0, 1)))) {
      markers <- markers[order(markers[, 1]), ]
    }
    chrNo <- unique(markers[, 1])
    for (i in 1:length(chrNo)) {
      markers[markers[, 1] == chrNo[i], ] <- markers[markers[,1] == chrNo[i], ][order(markers[markers[, 1] == chrNo[i], 2]), ]
    }
  }
  mkName <- rownames(markers)
  if (any(!is.element(rownames(genotypes), mkName))) {
    stop("Error: marker names do not match in markers and genotypes files")
  }
  if (any(rownames(genotypes) != mkName)) {
    loc2 <- NULL
    for (i in 1:nrow(genotypes)) {
      loc2 <- rbind(loc2, genotypes[which(rownames(genotypes) == mkName[i]), ])
    }
    rownames(loc2) <- mkName
    colnames(loc2) <- colnames(genotypes)
    genotypes <- loc2
  }
  indName <- colnames(genotypes)
  if (any(!is.element(colnames(traits), indName))) {
    stop("Error: individual names do not match in ", "genotypes and traits files")
  }
  if (any(colnames(traits) != indName)) {
    trait2 <- NULL
    for (i in 1:ncol(genotypes)) {
      trait2 <- rbind(trait2, traits[, which(colnames(traits) == indName[i])])
    }
    rownames(trait2) <- rownames(traits)
    colnames(trait2) <- indName
    traits <- trait2
  }
  if (!is.null(outputdir)) {
    dir <- dir.create(path = outputdir)
  }
  x <- y <- 1
  if (!is.null(outputdir)) {
    out <- c("output.txt", "qtlProfiles.csv"[is.null(qtlProfiles)], 
            "qtlSumm.csv"[is.null(qtlThres)], "corrZeroOrder.csv"[is.null(corrZeroOrder)], 
            "corrSecondOrder.csv"[is.null(corrSecondOrder)], 
            "corrPermutations.csv"[is.null(corrThres)], "network.sif"[cytoFiles], 
            "network.eda"[cytoFiles], "peakMultiplicity.csv"[!is.null(peaks)])
    fileExist <- list.files(path = outputdir)
    if (any(is.element(out, fileExist))) {
      cat("following files already exists in", outputdir, ": ", "\n\t")
      cat(out[which(is.element(out, fileExist))], sep = "\n\t")
      cat("Do you want to overwrite them: enter 1 for yes, 0 for no \n")
      x <- scan(n = 1)
    }
  }
  ck <- ls(envir = .GlobalEnv)
  ob <- c("qtlProfiles"[is.null(qtlProfiles)], "qtlSumm"[is.null(qtlThres)], 
        "fdr"[is.null(qtlThres)], "qtlThres"[is.null(qtlThres)], 
        "corrZeroOrder"[is.null(corrZeroOrder)], "corrSecondOrder"[is.null(corrSecondOrder)], 
        "corrPermutations"[is.null(corrThres)], "corrThres"[is.null(corrThres)])
  if (any(is.element(ob, ck))) {
    exist <- ob[which(is.element(ob, ck))]
    cat("Object(s) aleady exist: \n\t")
    cat(exist, sep = "\n\t")
    cat("Do you want to over-write them? enter 1 for yes, ")
    cat("enter 0 for no: \n")
    y <- scan(n = 1)
  }
  if (x != 1 | y != 1) {
    stop("Objects or files already exist. \n\t", "Please change the object names or specify other output directory")
  }else {
    if (is.vector(traits)) {
      traits <- t(as.matrix(traits))
    }
    name.traits <- rownames(traits)
    if (is.null(name.traits)) {
      name.traits <- 1:nrow(traits)
    }
    if (!is.null(outputdir)) {
      outfile <- paste(outputdir, "/", "output.txt", sep = "")
      cat("", file = outfile)
    }
    if (is.null(qtlProfiles)) {
      t1 <- proc.time()
      filename <- NULL
      cat("Step A: QTL mapping....", "\n")
      flush.console()
      if (!is.null(outputdir)) {
        filename <- paste(outputdir, "/qtlProfiles.csv", sep = "")
        cat("Step A: QTL mapping....", "\n", file = outfile, append = TRUE)
      }
      qtlProfiles <- qtlMapTwoPart(genotypes = genotypes, traits = traits, spike = spike, filename = filename)
      qtlProfiles <<- qtlProfiles
      t2 <- proc.time()
      cat("\t", "result in R object 'qtlProfiles'\n")
      cat("\t", "result in", filename, "\n")
      if (!is.null(outputdir)) {
        cat("\t", "result in R object 'qtlProfiles'\n", file = outfile, append = TRUE) 
        cat("\t", "result in", filename, "\n", file = outfile, append = TRUE)
      }
      timeToOutfile(t1,t2,outputdir,outfile)
    }else{
      cat("Step A: QTL mapping....skipped \n\t using user-provided QTL profiles\n\n")
      flush.console()
      if (!is.null(outputdir)) {
        cat("Step A: QTL mapping....skipped \n\t using user-provided QTL profiles", "\n", file = outfile, append = TRUE)
      }
    }
    
    if (is.null(qtlThres)) {
      n.simulations = 1000
      cat("Step B: Simulation test ( n = ", n.simulations, ") for QTL significance (-log10P) threshold....\n")
      flush.console()
      if (!is.null(outputdir)) {
        cat("Step B: Simulation test ( n =", n.simulations, ") for QTL significance (-log10P) threshold....\n", file = outfile, append = TRUE)
      }
      t1 <- proc.time()
      thres1 <- qtlThreshold(genotypes = genotypes, traits = traits, spike = spike, n.simulations = n.simulations, alpha = 0.05)
      fdr <- qtlFDR(qtlProfiles = qtlProfiles, fdrThres = 0.05, qtlThres = thres1)
      if (fdr[2, 1] > fdr[1, 1]) {
        thres1 <- fdr[1, 3]
        fdr2 <- fdr[1, 1]
      }else{
        fdr2 <- fdr[2, 1]
      }
      qtlThres <- round(thres1, digits = 2)
      qtlThres <<- qtlThres
      cat("\t", "alpha=0.05: QTL threshold = ", fdr[2,3], "\n")
      cat("\t", "fdr=0.05: QTL threshold = ", fdr[1, 3], "\n")
      cat("\t", "chose most stringent QTL threshold in R object 'qtlThres': \n")
      cat("\t", "logp =", qtlThres, "; FDR = ", fdr2, "\n")
      if (!is.null(outputdir)) {
        cat("\t", "alpha=0.05: qtl threshold = ", fdr[2,3], "\n", file = outfile, append = TRUE)
        cat("\t", "fdr=0.05: qtl threshold = ", fdr[1,3], "\n", file = outfile, append = TRUE)
        cat("\t", "chose most strigent threshold in R object 'qtlThres': \n", file = outfile, append = TRUE)
        cat("\t", "logp =", qtlThres, "; FDR = ", fdr2, "\n", file = outfile, append = TRUE)
      }
      t2 <- proc.time()
      timeToOutfile(t1,t2,outputdir,outfile)
    }else{
      cat("Step B: Simulation test for QTL significance threshold....skipped\n")
      cat("\t", "using user-provided QTL threshold:", qtlThres,"\n\n")
      flush.console()
      if (!is.null(outputdir)) {
        cat("Step B: Simulation test for QTL significance threshold....skipped\n", file = outfile, append = TRUE)
        cat("\t", "using user-provided QTL threshold:", qtlThres,"\n\n", file = outfile, append = TRUE)
      }
    }
    
    if(is.null(qtlSumm)){
      t1 <- proc.time()
      filename <- NULL
      cat("Step C: QTL summary....", "\n")
      flush.console()
      if (!is.null(outputdir)) {
        cat("Step C: QTL summary....", "\n", file = outfile, append = TRUE)
        filename <- paste(outputdir, "/qtlSumm.csv", sep = "")
      }
      qtlSumm <- qtlSummary(markers = markers, genotypes = genotypes, traits = traits, spike = spike, 
                 qtlProfiles = qtlProfiles, qtlThres = qtlThres, interval.dropoff = 1.5, filename = filename)
      qtlSumm <<- qtlSumm
      t2 <- proc.time()
      cat("\t", "result in R object: 'qtlSumm'", "\n")
      cat("\t", "result in", filename, "\n")
      if (!is.null(outputdir)) {
        cat("\t", "result in R object: 'qtlSumm'", "\n", file = outfile, append = TRUE)
        cat("\t", "result in", filename, "\n", file = outfile, append = TRUE)
      }
      timeToOutfile(t1,t2,outputdir,outfile)
    }else{
      cat("Step C: QTL summary....skipped \n\t using user-provided QTL summary\n\n")
      flush.console()
      if (!is.null(outputdir)) {
        cat("Step C: QTL summary....skipped \n\t using user-provided QTL summary\n\n", file = outfile, append = TRUE)
      }
    }
    
    if (is.null(corrZeroOrder)) {
      t1 <- proc.time()
      filename <- NULL
      cat("Step D: Zero-order correlation ....", "\n")
      flush.console()
      if (!is.null(outputdir)) {
        filename <- paste(outputdir, "/corrZeroOrder.csv", sep = "")
        cat("Step D: Zero-order correlation ....\n",file = outfile, append = TRUE)
      }
      if (corrMethod == "qtl") {
        corrZeroOrder <- qtlCorrZeroOrder(markers = markers, abs(qtlProfiles), qtlThres = qtlThres, filename = filename)
      }else if (corrMethod == "abundance") {
        corrZeroOrder <- cor(t(traits), method = "spearman", use = "pairwise.complete.obs")
      }
      corrZeroOrder <<- corrZeroOrder
      t2 <- proc.time()
      cat("\t", "result in R object: 'corrZeroOrder'\n")
      cat("\t", "result in", filename, "\n")
      if (!is.null(outputdir)) {
        cat("\t", "result in R object: 'corrZeroOrder'\n", file = outfile, append = TRUE)
        cat("\t", "result in", filename, "\n", file = outfile, append = TRUE)
      }
      timeToOutfile(t1,t2,outputdir,outfile)
    }else{
      cat("Step D: Zero-order correlation ....skipped \n\t using user-provided Zero-order correlation\n\n")
      flush.console()
      if (!is.null(outputdir)) {
        cat("Step D: Zero-order correlation ....skipped \n\t using user-provided Zero-order correlation\n\n", file = outfile, append = TRUE)
      }
    }
    
    if (is.null(corrSecondOrder)) {
      t1 <- proc.time()
      filename <- NULL
      cat("Step E: 2nd-order correlation ....", "\n")
      flush.console()
      if (!is.null(outputdir)) {
        filename <- paste(outputdir, "/corrSecondOrder.csv", sep = "")
        cat("Step E: 2nd-order correlation of QTL profiles....\n", file = outfile, append = TRUE)
      }
      corrSecondOrder <- qtlCorrSecondOrder(corrZeroOrder, filename = filename)
      corrSecondOrder <<- corrSecondOrder
      t2 <- proc.time()
      cat("\t", "result in R object: 'corrSecondOrder'\n")
      cat("\t", "result in", filename, "\n")
      if (!is.null(outputdir)) {
          cat("\t", "result in R object: 'corrSecondOrder'\n", file = outfile, append = TRUE)
          cat("\t", "result in", filename, "\n", file = outfile, append = TRUE)
      }
      timeToOutfile(t1,t2,outputdir,outfile)
    }else{
      cat("Step E: 2nd-order correlation ....skipped \n\t using user-provided 2nd-order correlation\n\n")
      flush.console()
      if (!is.null(outputdir)) {
        cat("Step E: 2nd-order correlation ....skipped \n\t using user-provided 2nd-order correlation\n\n", file = outfile, append = TRUE)
      }
    }
    
    if (is.null(corrThres)) {
      n.permutations <- 10000
      t1 <- proc.time()
      cat("Step F: Permutation test ( n.permutations = ", n.permutations, ") for 2nd-order correlation significance threshold...\n")
      flush.console()
      if (!is.null(outputdir)) {
        cat("Step F: permutation test ( n.permutations = ", n.permutations, ") for 2nd-order correlation significance threshold...\n", file = outfile, append = TRUE)
      }
      corrThres <- qtlCorrThreshold(markers = markers, genotypes = genotypes, traits = traits, spike = spike, qtlThres = qtlThres, method = corrMethod, n.permutations = n.permutations)
      corrThres <<- corrThres
      t2 <- proc.time()
      cat("\t", "result in R object: 'corrThres'", "\n")
      cat("\t", "value of each permutation in R object: 'corrPermutations'\n")
      cat("\t", "Chosen correlation threshold (alpha = 0.05):",corrThres, "\n")
      if (!is.null(outputdir)) {
        cat("\t", "result in R object: 'corrThres'\n", file = outfile, append = TRUE)
        cat("\t", "Chosen correlation threshold (alpha = 0.05):",corrThres, "\n", file = outfile, append = TRUE)
      }
      timeToOutfile(t1,t2,outputdir,outfile)
    }else{
      cat("Step F: Permutation test for 2nd-order correlation significance threshold...skipped\n")
      cat("\t", "using user-provided correlation threshold:", corrThres, "\n\n")
      flush.console()
      if (!is.null(outputdir)) {
        cat("Step F: Permutation test for 2nd-order correlation significance threshold...skipped\n", file = outfile, append = TRUE)
        cat("\t", "using user-provided correlation threshold:", corrThres, "\n\n", file = outfile, append = TRUE)
      }
    }
      
    if (cytoFiles & !is.null(outputdir)) {
      cat("Step G: Create Cytoscape network files...\n")
      filename <- paste(outputdir, "/network", sep = "")
      createCytoFiles(corrSecondOrder, filename = filename, simThres = corrThres, hideNodes = TRUE)
      cat("\t", "SIF file is:", paste(filename, ".sif", sep = ""), "\n")
      cat("\t", "EDA file is:", paste(filename, ".eda", sep = ""), "\n\n")
      cat("Step G: Create Cytoscape network files...\n", file = outfile, append = TRUE)
      cat("\t", "SIF file is:", paste(filename, ".sif", sep = ""), "\n", file = outfile, append = TRUE)
      cat("\t", "EDA file is:", paste(filename, ".eda", sep = ""), "\n\n", file = outfile, append = TRUE)
    }else{
      cat("Step G: Create Cytoscape network files...skipped \n\t because cytoFiles = FALSE or outputdir = NULL\n\n")
      if (!is.null(outputdir)) {
        cat("Step G: Create Cytoscape network files...skipped because cytoFiles = FALSE\n\n", file = outfile, append = TRUE)
      }
    }
      
    if (!is.null(peaks)) {
      t1 <- proc.time()
      filename <- NULL
      cat("Step H: Find peak multiplicity ....", "\n")
      flush.console()
      if (!is.null(outputdir)) {
          cat("Step H: Find peak multiplicity ....", "\n", file = outfile, append = TRUE)
          filename <- paste(outputdir, "/peakMultiplicity.csv", sep = "")
      }
      peakMultiplicity <- findPeakMultiplicity(corrZeroOrder, peaks, corrThres = 0.95, filename = filename)
      peakMultiplicity <<- peakMultiplicity
      t2 <- proc.time()
      cat("\t", "result in R object: 'peakMultiplicity'\n")
      cat("\t", "result in", filename, "\n")
      if (!is.null(outputdir)) {
          cat("\t", "result in R object: 'peakMultiplicity'\n", file = outfile, append = TRUE)
          cat("\t", "result in", filename, "\n", file = outfile, append = TRUE)
      }
      timeToOutfile(t1,t2,outputdir,outfile)
    }else{
        cat("Step H: Detection of peak multiplicity...skipped\n\n")
        flush.console()
        if (!is.null(outputdir)) {
            cat("Step H: Detection of peak multiplicity...skipped\n", file = outfile, append = TRUE)
        }
    }
  }
}