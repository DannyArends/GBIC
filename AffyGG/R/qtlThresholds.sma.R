# Function name: qtlThresholds.sma
# Author: Rudi Alberts <r.alberts@rug.nl>
# Maintainer: Gonzalo Vera <g.vera.rodriguez@rug.nl>
# Version: 1.2.0
# Date: 1 Feb. 2006

qtlThresholds.sma <- function( genotypes, batch, 
                               nProbes=16, nIndiv=30, 
                               n.simulations=10000, filename=NULL )
{
  library(nlme)                 # necessary for lme function

  nValues <- nProbes * nIndiv   # number of values necessary per probe set
  
  # generate one vector with random expression data with mean 0 and sd 1.
  rand    <- rnorm( nValues*n.simulations )                                    

  # We preset a general part for anova input. Marker will be added
  anovainalg <- data.frame(     # Three columns: batch, probe and mouse
              batch= rep( batch, each=nProbes ),
              probe= 1:nProbes, 
              mouse= rep( 1:nIndiv  , each=nProbes ) ) 

  # vectors to store simulation results:
  maxFmarker      <- NULL       # F-value QTL
  maxFinteraction <- NULL       # F-value QTL:probe
  maxpmarker      <- NULL       # P-value QTL
  maxpinteraction <- NULL       # P-value QTL:probe

  for( idxSimul in 0:(n.simulations-1) ) # simulate 10000 times
  {
    # get data for one simulation from the rand vector
    exprR <- round( rand[ (idxSimul*nValues+1):(idxSimul*nValues+nValues) ], 2 )
    Fmarker       <- NULL # vectors to store results of this simulation
    Finteraction  <- NULL
    pmarker       <- NULL
    pinteraction  <- NULL

    for( idxMarker in 1:nrow(genotypes) )  # loop over all markers
    {
      # cat(mar)
      
      # add marker and expression data to anovain
      anovain <- data.frame( anovainalg, 
                             marker=rep( as.numeric(genotypes[idxMarker,]), each=nProbes),
                             value=exprR )
      # fit linear mixed-effects model
      g <- lme( value~factor(batch) + factor(probe) +
                factor(batch):factor(probe) + factor(marker) +
                factor(marker):factor(probe),
                data=anovain,
                random=~1|factor(mouse) )
      ag <- anova(g) # compute anova table
      # extract F and p values from anova table
      Fmarker       <- c( Fmarker,      as.numeric(ag[[3]][4]))
      Finteraction  <- c( Finteraction, as.numeric(ag[[3]][6]))
      pmarker       <- c( pmarker,      as.numeric(ag[[4]][4]))
      pinteraction  <- c( pinteraction, as.numeric(ag[[4]][6]))
    }
    # collect maximum F and p values for all simulations
    maxFmarker      <-c( maxFmarker,      max(Fmarker) )
    maxFinteraction <-c( maxFinteraction, max(Finteraction) )
    maxpmarker      <-c( maxpmarker,      max(pmarker) )
    maxpinteraction <-c( maxpinteraction, max(pinteraction) )
  }
  # collect all results in one data.frame
  thresholds <- data.frame( maxFmarker=maxFmarker, maxFinteraction=maxFinteraction,
                            maxpmarker=maxpmarker, maxpinteraction=maxpinteraction )
  
  if (!is.null(filename)) 
  {
     write.table( thresholds, file=filename, sep="," ) # store results
  }
  
  thresholds
}