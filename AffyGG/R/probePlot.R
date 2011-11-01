# Function name: probePlot
# Author: Rudi Alberts <r.alberts@rug.nl>
# Maintainer: Gonzalo Vera <g.vera.rodriguez@rug.nl>
# Version: 1.1.0
# Date: 1 Feb. 2006

probePlot <- function( traits, probesetName, probesPos, alleleColors=NULL, 
                       filename=NULL ) 
{
  if( !is.null( filename ) )
  { 
    png( filename, height = 550, width = 700 ) # open png file
  }
  if( is.null( alleleColors ) )
  {
   alleleColors <- 1                           # Black color assigned by default
  }
  
  
  probesPos <- probesPos - min(probesPos) + 1  # let positions start from 1
  lengthprobe <- 25
  par( mai=c(0, 1, 0.2, 0.3) )                 # set margins
  par( fig=c(0,1,0.84,1) )                     # upper part of figure
  
  # plot empty box
  plot( c( probesPos, max(probesPos)+lengthprobe+2 ), 
        c( 1:length(probesPos), length(probesPos) ),
        xlab= paste( "probeset", probesetName ),
        ylab= "probe", type="n", xaxt="n" )
  
  # per probe, plot line and probe number
  for( i in 1:length(probesPos) ) 
  {
    lines( c(probesPos[i], probesPos[i]+lengthprobe), c(i,i), lwd=1, col="blue" )
    text( probesPos[i]+13, i+0.2, i )
  }
  
  text( probesPos[length(probesPos)]+25, 2, "Relative probe position (bp)", pos=2)
  par( mai=c(1.2, 1, 0.1, 0.3) )              # set margins
  par( fig=c(0,1,0,0.8), new=T, xpd=T )       # lower part of figure
  myxlab <- paste( "probe signals for probeset", probesetName )   
  
  # plot the signals (either PM or MM)
  matplot( traits, type="l", main="", xlab=myxlab, xaxt="n", 
           ylab="log2 intensity", lty=1:ncol( traits ), ylim=range(5:15),
           col=alleleColors)
  axis( 1, 1:nrow(traits), 1:nrow(traits) ) # add numbers to X axis
  
  if( !is.null( filename ) )
  {
    dev.off()
  }
}

