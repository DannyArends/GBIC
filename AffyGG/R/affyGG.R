# Function name: affyGG
# Author: Rudi Alberts <r.alberts@rug.nl>
# Maintainer: Gonzalo Vera <g.vera.rodriguez@rug.nl>
# GIT Maintainer: Danny Arends <d.arends@gmail.com>
# Version: 1.0.1
# Date: 1 May. 2007

affyGG <- function( genotypes, batch, nProbes=16, nIndiv=30, n.simulations=100,
                    directory=NULL, celfiles, probelevel=T,
                    probesetName=NULL, overwrite=F, n.Plots=5,
                    markersPos=NULL, probesPos=NULL, chrOffsets=NULL ){

  # 0. Welcome
  cat( "\nPackage: affyGG \nComputational protocols for genetical genomics analyses with Affymetrix arrays. \n\n" )
  initTime <- proc.time()

  # 1. Checking arguments... 
  cat( "STEP 1: Checking arguments... " )  
  flush.console()

  # The Working Directory is setup at the same path as the .CEL files
  if( is.null( directory ) ) {
    stop( "A working directory with .CEL files should be provided")
  }
  else {
    setwd( directory )
  }
  # We need to make that variable, function scope' accesible
  qtlMap <- NULL
  
  # Non-optional arguments should be provided
  if( is.null( markersPos ) ) {
    stop( "A value for argument: markerPos should be provided")
  }
  if( is.null( chrOffsets ) ) {
    stop( "A value for argument: chrOffsets should be provided")
  }
  if( is.null( probesPos ) ) {
    stop( "A value for argument: probePos should be provided")
  }
  
  # We show time elapsed
  prevTime <- partialTime <- proc.time()
  cat( " (T: ", (partialTime - initTime)[3], ", P: ", (prevTime - initTime)[3], 
       " secs elapsed) \n", sep="" )
  flush.console()
  
  # 2. Preprocessing Affymetrix .CEL files
  cat( "STEP 2: Preprocessing Affymetrix .CEL files... " )
  flush.console()
  
  #   We check if the output file exists, so, it was previously calculated
  if( file.exists("probesignals.csv") && !overwrite )
  {
    # We load the data previously calculated
    probesignals <- read.csv( "probesignals.csv" )
  }
  #   If not, we had to process the .CEL files again
  else 
  {
    probesignals <- rma.preprocessing( directory, celfiles, "probesignals.csv" )
  }
  # We show time elapsed
  partialTime <- proc.time()
  cat( " (T: ", (partialTime - initTime)[3], ", P: ", (partialTime - prevTime)[3], 
       " secs elapsed) \n", sep="" )
  prevTime <- partialTime 
  flush.console()

  
  # 3.1 Significance thresholds calculation
  cat( "STEP 3.1: Significance thresholds calculation... " )
  flush.console()    
  
  #   We check if the output file exists, so, it was previously calculated    
  if( file.exists("qtlthresholds.csv") && !overwrite )
  {
    qtlThres_temp <- read.csv("qtlthresholds.csv")
  }
  else
  {
    qtlThres_temp <- qtlThresholds.sma( genotypes, batch, nProbes, nIndiv,
                                        n.simulations, "qtlthresholds.csv")
  }
  # We show time elapsed
  partialTime <- proc.time()
  cat( " (T: ", (partialTime - initTime)[3], ", P: ", (partialTime - prevTime)[3], 
       " secs elapsed) \n", sep="" )
  prevTime <- partialTime 
  flush.console()

  # 3.2 Choosing the 95 percentile value
  cat( "STEP 3.2: Choosing the 95 percentile value... " )
  flush.console()
  
  qtlThres <- -log10(sort( qtlThres_temp$maxpmarker, 
                           decreasing=T )[round(n.simulations*.95)])
  intThres <- -log10(sort( qtlThres_temp$maxpinteraction, 
                           decreasing=T )[round(n.simulations*.95)])
  # We show time elapsed
  partialTime <- proc.time()
  cat( " (T: ", (partialTime - initTime)[3], ", P: ", (partialTime - prevTime)[3], 
       " secs elapsed) \n", sep="" )
  prevTime <- partialTime 
  flush.console()

  # 4.1 QTL analysis
  cat( "STEP 4.1: QTL analysis... " ) 
  flush.console()
  
  listProbeSets    <- as.character(unique( probesignals$probeset ))
  if( is.null( probesetName ) ) 
  {
     maxQtl<-NULL
     for( idxProbeSet in 1:length( listProbeSets ) )
     {
        traits <- probesignals[ probesignals$probeset == listProbeSets[ idxProbeSet ] , ]
        traits <- traits[,-1]  # strip the probe set name from traits
        nameProbeSet <- auxCleanFilename( listProbeSets[ idxProbeSet ] )
        if( probelevel )  {
          qtlMapFile   <- paste( "qtlMap.xProbe.", nameProbeSet, ".csv", sep="" )
          qtlMap <- qtlMap.xProbe( genotypes, traits, batch, qtlMapFile )
          maxQtl <- rbind( maxQtl, data.frame(probeset = listProbeSets[ idxProbeSet ],
                           maxQtl = min(qtlMap$pmarkerperprobe) ,
                           markerName = qtlMap$marker[which.min(qtlMap$pmarkerperprobe)]))
        }
        else {
          qtlMapFile   <- paste( "qtlMap.xProbeSet.", nameProbeSet, ".csv", sep="" )
          qtlMap <- qtlMap.xProbeSet( genotypes, traits, batch, qtlMapFile )
          maxQtl <- rbind( maxQtl, data.frame(probeset = listProbeSets[ idxProbeSet ],
                           maxQtl = min(qtlMap$pmarker),
                           markerName = qtlMap$markers[which.min(qtlMap$pmarker)]))
        }
     }
  }
  else   # probesetName NOT NULL
  {
    # We check if the probesetName provided, exists.
    if( probesetName %in% listProbeSets )
    {
      traits <- probesignals[ probesignals$probeset == probesetName , ]
      traits <- traits[,-1]  # strip the probe set name from traits
      nameProbeSet <- auxCleanFilename( probesetName )

      if( probelevel )  {
          qtlMapFile   <- paste( "qtlMap.xProbe.", nameProbeSet, ".csv", sep="" )
          qtlMap <- qtlMap.xProbe( genotypes, traits, batch, qtlMapFile )
          maxQtl <- data.frame(probeset = probesetName,
                           maxQtl = min(qtlMap$pmarkerperprobe) ,
                           markerName = qtlMap$marker[which.min(qtlMap$pmarkerperprobe)])
      }
      else {
        qtlMapFile   <- paste( "qtlMap.xProbeSet.", nameProbeSet, ".csv", sep="" )
        qtlMap <- qtlMap.xProbeSet( genotypes, traits, batch, qtlMapFile )
        maxQtl <- data.frame(probeset = probesetName,
                           maxQtl = min(qtlMap$pmarker),
                           markerName = qtlMap$markers[which.min(qtlMap$pmarker)])
      }
    }
    else {
      stop("The probeset name indicated couldn't be found.")
    }
  }
  # We show time elapsed
  partialTime <- proc.time()
  cat( " (T: ", (partialTime - initTime)[3], ", P: ", (partialTime - prevTime)[3], 
       " secs elapsed) \n", sep="" )
  prevTime <- partialTime 
  flush.console()
  
  # 4.2. Sort QTL results  
  if( is.null( probesetName ) ) 
  {
    cat( "STEP 4.2: Sorting QTL results... " )
    flush.console()
  
    maxQtl <- maxQtl[order(maxQtl$maxQtl),]
    # We show time elapsed
    partialTime <- proc.time()
    cat( " (T: ", (partialTime - initTime)[3], ", P: ", (partialTime - prevTime)[3], 
         " secs elapsed) \n", sep="" )
    prevTime <- partialTime 
    flush.console()
  }
  
  # 5. Plot probe signals
  cat( "STEP 5: Plotting probe signals... ")
  flush.console()
  
  # 5.1 If the user doesn't provide a probesetName, we plot the first n.Plots
  if( is.null( probesetName ) )
  {
    for (idxProbeSet in 1:n.Plots)
    {
          nameProbeSet <- auxCleanFilename( maxQtl$probeset[idxProbeSet] )
          traits <- probesignals[ probesignals$probeset == maxQtl$probeset[idxProbeSet] , ]
          alleleColors <- as.numeric( genotypes[maxQtl$markerName[idxProbeSet],] )
          probePlotFile   <- paste( "probePlot.", nameProbeSet, ".png", sep="" )
          probePlot( traits, maxQtl$probeset[idxProbeSet], probesPos, alleleColors, probePlotFile )
    }
  }
  # 5.2 If probesetName is provided, then we only plot once
  else
  {
    nameProbeSet <- auxCleanFilename( probesetName )
    traits <- probesignals[ probesignals$probeset == probesetName , ]
    alleleColors <- as.numeric( genotypes[maxQtl$markerName,] )
    probePlotFile   <- paste( "probePlot.", nameProbeSet, ".png", sep="" )
    probePlot( traits, probesetName, probesPos, alleleColors, probePlotFile )
  }
  # We show time elapsed
  partialTime <- proc.time()
  cat( " (T: ", (partialTime - initTime)[3], ", P: ", (partialTime - prevTime)[3], 
       " secs elapsed) \n", sep="" )
  prevTime <- partialTime 
  flush.console()


  # 6. Make QTL plots
  cat( "STEP 6: Plotting QTL maps... ")
  flush.console()
  
  # 6.1 If the user doesn't provide a probesetName, we plot the first n.Plots
  if( is.null( probesetName ) )
  {
    for (idxProbeSet in 1:n.Plots)
    {
      nameProbeSet <- auxCleanFilename( maxQtl$probeset[idxProbeSet] )
      if ( probelevel )
      {
        probeQtlProfiles <- read.csv( paste( "qtlMap.xProbe.",
                                              nameProbeSet, ".csv", sep="" ) )
        qtlPlotFile   <- paste( "qtlPlot.xProbe.", nameProbeSet, ".png", sep="" )
        qtlPlot.xProbe( probesetName, markersPos, probeQtlProfiles,
                           qtlThres, chrOffsets, qtlPlotFile )
      }
      else {

        bothProfiles <- read.csv( paste( "qtlMap.xProbeSet.", nameProbeSet, ".csv", sep="" ) )
        probesetQtlProfile <- -log10( bothProfiles$pmarker )
        interactionProfile <- -log10( bothProfiles$pinteraction )
        interactionThres <- intThres
        
        qtlPlotFile   <- paste( "qtlPlot.xProbeSet.", nameProbeSet, ".png", sep="" )
        qtlPlot.xProbeSet( probesetName, markersPos,
                      probesetQtlProfile, interactionProfile,
                      qtlThres, interactionThres,
                      chrOffsets, qtlPlotFile )
      }
    }
  }
  #6.2 If probesetName is provided, then we only plot once
  else
  {   
    nameProbeSet <- auxCleanFilename( probesetName )
    if ( probelevel )
    {
      probeQtlProfiles <- qtlMap
      qtlPlotFile   <- paste( "qtlPlot.xProbe.", nameProbeSet, ".png", sep="" )
      qtlPlot.xProbe( probesetName, markersPos, probeQtlProfiles,
                         qtlThres, chrOffsets, qtlPlotFile )
    }
    else 
    {   
      bothProfiles <- qtlMap
      probesetQtlProfile <- -log10( bothProfiles$pmarker )
      interactionProfile <- -log10( bothProfiles$pinteraction )
      interactionThres <- intThres

      qtlPlotFile   <- paste( "qtlPlot.xProbeSet.", nameProbeSet, ".png", sep="" )
      qtlPlot.xProbeSet( probesetName, markersPos,
                    probesetQtlProfile, interactionProfile,
                    qtlThres, interactionThres,
                    chrOffsets, qtlPlotFile )
    }
  }
  # We show time elapsed
  partialTime <- proc.time()
  cat( " (T: ", (partialTime - initTime)[3], ", P: ", (partialTime - prevTime)[3], 
       " secs elapsed) \n", sep="" )
  prevTime <- partialTime 
  flush.console()       
}
