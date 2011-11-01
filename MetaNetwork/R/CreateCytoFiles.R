#
# CreateCytoFiles.R
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

createCytoFiles <- function (corrMatrix, filename, simThres = NULL, hideNodes = TRUE){
    disThres <- NULL
    relationship <- "cc"
    name <- rownames(corrMatrix)
    if (is.null(name)) {
        name <- paste("node", 1:nrow(corrMatrix), sep = "")
    }
    sif.file <- paste(filename, ".sif", sep = "")
    eda.file <- paste(filename, ".eda", sep = "")
    if (is.null(simThres) & is.null(disThres)) {
        stop("no similarityThreshold for network edges")
    }
    cat("", file = sif.file)
    cat("InteractionStrength", "\n", sep = "", file = eda.file)
    for (i in 2:nrow(corrMatrix)) {
        for (j in 1:(i - 1)) {
            if (!is.null(simThres)) {
                if (abs(corrMatrix[i, j]) > simThres) {
                  cat(name[i], relationship, name[j], "\n", sep = " ", file = sif.file, append = TRUE)
                  cat(name[i], paste("(", relationship, ")", sep = ""), name[j], "=", corrMatrix[i, j], "\n", sep = " ", file = eda.file, append = TRUE)
                }else if (!hideNodes) {
                  cat(name[i], "\n", sep = " ", file = sif.file, append = TRUE)
                }
            }
            else if (!is.null(disThres)) {
                if (abs(corrMatrix[i, j]) < disThres) {
                  cat(name[i], relationship, name[j], "\n", sep = " ", file = sif.file, append = TRUE)
                  cat(name[i], paste("(", relationship, ")", sep = ""), name[j], "=", corrMatrix[i, j], " ", sep = "\t", file = eda.file, append = TRUE)
                }
                else if (!hideNodes) {
                  cat(name[i], "\n", sep = " ", file = sif.file, append = TRUE)
                }
            }
        }
    }
}