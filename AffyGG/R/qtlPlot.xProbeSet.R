# Function name: qtlPlot.xProbeSet
# Author: Rudi Alberts <r.alberts@rug.nl>
# Maintainer: Gonzalo Vera <g.vera.rodriguez@rug.nl>
# Version: 1.0.0
# Date: 1 Feb. 2006

qtlPlot.xProbeSet <- function( probesetName, markersPos, 
                               probesetQtlProfile, interactionProfile,
                               qtlThres, interactionThres, 
                               chrOffsets, filename=NULL ) 
{
  if( !is.null( filename ) )
  {
    png( filename=filename, bg="white", width=1100, height=300 )
  }

  markersPos <- markersPos[,2]  # We only use the bp column
  nChr <- length( chrOffsets ) 
  par( mai=c(1, 1, 0.5, 0.3) )                                # set margins  
  maxheight <- max( probesetQtlProfile, interactionProfile )  # get max for scaling the plot
  
  # make empty plot of wished size
  plot( markersPos, probesetQtlProfile, 
        type= "n", 
        main= paste("QTL analysis per probe set for probe set", probesetName),
        xlab= "marker", 
        ylab= "-10log(p value)", 
        ylim= range(0:maxheight*1.2))  
  
  for( i in 1:nChr ) 
  { 
    # add lines and names for chr's
    lines( c( chrOffsets[i], chrOffsets[i] ), 
           c( 0,maxheight ), lty=1, lwd=1, col="grey") 
    
    if( i < nChr  ) 
    { 
      text( chrOffsets[i], maxheight, i, cex=1.2, pos=4 ) 
    }
    else if( i == nChr ) 
    { 
      lines( c(max(markersPos),max(markersPos)), 
             c(0,maxheight), lty=1, lwd=1, col="grey") 
      text( chrOffsets[i], maxheight, "X", cex=1.2, pos=4)
    }
  }
  
  # plot -10log(p) values and thresholds
  lines( markersPos, interactionProfile, col="green4") 
  lines( markersPos, probesetQtlProfile, col="blue")
  lines( c(0,max(markersPos)), 
         c(interactionThres, interactionThres), col="green4", lty=2)
  lines( c(0,max(markersPos)), c(qtlThres, qtlThres), col="blue", lty=2)
  
  for (idxMarkers in 1:length(markersPos)) 
  {
     lines( c( markersPos[idxMarkers],markersPos[idxMarkers] ),
            c( -maxheight/100,-maxheight/25 ), col="darkgrey" )
  }
  
  if( !is.null( filename ) )
  {
    dev.off()
  }
}

