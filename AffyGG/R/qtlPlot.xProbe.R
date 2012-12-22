# Package: GeneticalGenomics v1.0.0
# Originally part of Rudi Alberts' PhD Thesis
# Function name: gg.qtlPlot.xProbe
# Author: Rudi Alberts <r.alberts@rug.nl>
# Maintainer: Gonzalo Vera <g.vera.rodriguez@rug.nl>
# Version: 1.0.0
# Date: 1 Feb. 2006

qtlPlot.xProbe <- function( probesetName, markersPos, probeQtlProfiles,
                               qtlThres, chrOffsets, filename=NULL ) 
{
  if( !is.null( filename ) )
  {
    png(filename=filename, bg="white", width=1100, height=700)
  }
  
  markersPos <- markersPos[,2]  # We only use the bp column
  nChr <- length( chrOffsets ) 
  nProbes <- length(unique(probeQtlProfiles$probenr)) # get number of probes
  maxheight <- max( probeQtlProfiles[,4] ) # get max for scaling the plot
  stepsize  <- maxheight
  
  # make empty plot of wished size
  plot( markersPos, probeQtlProfiles[ probeQtlProfiles$probenr==1, 4], 
        type="n", yaxt="n",
        main=paste("QTL analysis per probe for probe set",probesetName),
        xlab="marker", ylab="-10log(p)", ylim=range(0:(maxheight*110*nProbes)/100),
        cex=1) # make empty plot of wished size
  
  totmaxheight <- (maxheight*110*nProbes)/100

  for (i in 1:nChr) 
  {
    lines( c(chrOffsets[i], chrOffsets[i]), c(0,totmaxheight), 
           lty=1, lwd=1, col="grey" ) # add lines and names for chr's
    if( i==1 ) {
        text(chrOffsets[i], totmaxheight,"chr 1", pos=4)
    }
    else if( i < nChr && i > 1 ) {
        text(chrOffsets[i], totmaxheight, i, pos=4)
    }
    else if( i == nChr ) {
      lines( c(max(markersPos),max(markersPos)), 
             c(0,totmaxheight), lty=1, lwd=1, col="grey" )
      text( chrOffsets[i], totmaxheight, "X", pos=4 )      
    }
  }
  
  mycols <- rep( c("orange1","orange4","lightblue4","blue4"), 6 )  
  for( i in 1:nProbes ) 
  {
    # plot -10log(p) values and thresholds per probe
    pp <- data.frame( probeQtlProfiles[probeQtlProfiles$probenr==i, 4] ) 
    lines( markersPos, pp[,1]+(i-1)*stepsize, col=mycols[i] )  
    text( 24, (i-1)*stepsize, i, col=mycols[i], pos=2 )
    lines( c( 0, max(markersPos) ), 
           c( (i-1)*stepsize + qtlThres, (i-1)*stepsize + qtlThres),
          lty=2, col=mycols[i] )
  }
  for (idxMarkers in 1:length(markersPos)) 
  {
     lines( c( markersPos[idxMarkers],markersPos[idxMarkers] ),
            c( .4-stepsize,5-stepsize ), col="darkgrey" )
  }
  
  
  if( !is.null( filename ) )
  {
    dev.off()
  }
}
