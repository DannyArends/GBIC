# Function name: splitProbeSetTraits
# Author: Gonzalo Vera <g.vera.rodriguez@rug.nl>
# Maintainer: Gonzalo Vera <g.vera.rodriguez@rug.nl>
# Version: 1.0.0
# Date: 1 Feb. 2006

splitProbeSetTraits <- function(affy.signals, directory){
  listProbeSets    <- unique(affy.signals$probeset)

  nameIndividuals  <- toupper(colnames( affy.signals )[-1])
  nameIndividuals  <- strsplit(nameIndividuals, ".CEL")


  for( idxProbeSet in 1:length( listProbeSets )){
    # Check that future file name doesn't contains reserved caracthers like ? / \ = < > : ; " ,
    nameProbeSet <- auxCleanFilename(listProbeSets[ idxProbeSet ])
    
    nProbes    <- nrow( affy.signals[ affy.signals$probeset ==  listProbeSets[idxProbeSet] , ] )
    traitsTemp <- data.frame( affy.signals[ affy.signals$probeset ==  listProbeSets[idxProbeSet] , ][,-1],
                                row.names= paste(nameProbeSet, 1:nProbes, sep=""))
    colnames(traitsTemp) <- nameIndividuals

    traitsFile <- paste( directory, "/traits.", nameProbeSet, ".csv", sep="" )
    write.csv( traitsTemp, file=traitsFile )
  }
}

