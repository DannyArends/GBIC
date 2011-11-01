# Function name: gg.qtlMap.xProbe
# Author: Rudi Alberts <r.alberts@rug.nl>
# Maintainer: Gonzalo Vera <g.vera.rodriguez@rug.nl>
# Version: 1.2.0
# Date: 1 Feb. 2006

qtlMap.xProbe <- function( genotypes, traits, batch=NULL, filename=NULL )
{
  nProbes <- nrow( traits )           # number of probes
  onevec  <- c( as.matrix( traits ) ) # put data in one vector
                           
  anovainalg <- data.frame( probe= 1:nProbes,
                            mouse= rep(1:ncol(traits), each=nProbes), 
                            value= onevec ) # general
                            # part of anova input, only marker will be added
  if (!is.null(batch)) 
  {
    anovainalg <- data.frame( batch= rep( batch, each=nProbes), anovainalg )  
  }

 # vector to store results:
  perprobe        <- NULL   # results

  for( idxMarker in 1:nrow( genotypes ) )  #  QTL mapping on all markers
  {
    Fmarkerperprobe <- NULL   # F-value QTL
    pmarkerperprobe <- NULL   # p-value QTL
    
    # add marker data to anovainnow
    anovainnow <- cbind( anovainalg, 
                         marker=rep(as.numeric(genotypes[idxMarker,]), each=nProbes) )
    for( idxProbe in 1:nProbes )
    {
      # select data for probe p
      nowin <- anovainnow[ anovainnow$probe == idxProbe,] 
      
      # fit linear model
      if (!is.null(batch)) {
        g <- lm( value~factor(batch)+factor(marker),nowin ) 
      }
      else {
        g <- lm( value~factor(marker),nowin )       
      }
      
      # compute anova table
      ag    <- anova(g) 
      
      # get F and p values from anova
      if (!is.null(batch)) {
        Fmarkerperprobe <- c( Fmarkerperprobe,as.numeric(ag[[4]][2]))
        pmarkerperprobe <- c( pmarkerperprobe,as.numeric(ag[[5]][2]))
      }
      else {      
        Fmarkerperprobe <- c( Fmarkerperprobe,as.numeric(ag[[4]][1]))
        pmarkerperprobe <- c( pmarkerperprobe,as.numeric(ag[[5]][1]))
      }
    }

    temp <- data.frame( marker=rownames(genotypes)[idxMarker], 
                        probenr=1:nProbes,
                        Fmarkerperprobe,
                        pmarkerperprobe )
    rownames(temp) <- paste( temp$marker, 1:nProbes, sep="" )
    perprobe       <- rbind( perprobe,temp ) # collect results
  }

  rownames(perprobe) <- 1:nrow(perprobe)
  
  if (!is.null(filename))    
  {
    write.table( perprobe, file=filename, sep="," )
  }
  
  return(perprobe) #return results
}
