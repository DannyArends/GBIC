# Function name: gg.qtlMap.xProbeSet
# Author: Rudi Alberts <r.alberts@rug.nl>
# Maintainer: Gonzalo Vera <g.vera.rodriguez@rug.nl>
# Version: 1.2.0
# Date: 1 Feb. 2006

qtlMap.xProbeSet <- function( genotypes, traits, batch=NULL, filename=NULL )
{
  library(nlme)                 # necessary for lme function

  nProbes <- nrow( traits )  # number of probes
  onevec <- c( as.matrix( traits ) ) # put data in one vector    

  anovainalg <- data.frame( probe= 1:nProbes,
                            mouse= rep(1:ncol( traits ), each=nProbes),
                            value= onevec ) # general part of anova input,
                            # only marker will be added
  if (!is.null(batch)) 
  {
    anovainalg <- data.frame( batch= rep( batch, each=nProbes), anovainalg )  
  }
  # vectors to store results:
  Fmarker       <- NULL      # F-value QTL
  Finteraction  <- NULL      # F-value QTL:probe
  pmarker       <- NULL      # P-value QTL
  pinteraction  <- NULL      # P-value QTL:probe

  for (idxMarker in 1:nrow(genotypes)) #  QTL mapping on all markers
  {
    # add marker data to anovainnow

    anovainnow <- cbind( anovainalg, 
                         marker=rep(as.numeric(genotypes[idxMarker,]), each=nProbes) )
    
    # fit linear mixed-effects model
    if (!is.null(batch)) {    
      g   <- lme( value~factor(batch) + factor(probe) +
                factor(batch):factor(probe) + factor(marker) +
                factor(marker):factor(probe),
                data=anovainnow,
                random=~1|factor(mouse) )
    }
    else {
       g   <- lme( value~factor(probe) + factor(marker) +
                factor(marker):factor(probe),
                data=anovainnow,
                random=~1|factor(mouse) )
    }
                
    # compute anova table
    ag  <- anova(g) 
    
    if (!is.null(batch)) {    
      # get F and p values from anova    
      Fmarker      <- c( Fmarker,      as.numeric(ag[[3]][4]))
      Finteraction <- c( Finteraction, as.numeric(ag[[3]][6]))
      
      # derive p values using the pf function, to avoid zero p values of which
      # you can not take a log
      pmarker      <- c( pmarker,      pf(1/as.numeric(ag[4,3]),ag[4,2],ag[4,1]) )
      pinteraction <- c( pinteraction, pf(1/as.numeric(ag[6,3]),ag[6,2],ag[6,1]) )
    }
    else {
      # get F and p values from anova    
      Fmarker      <- c( Fmarker,      as.numeric(ag[[3]][3]))
      Finteraction <- c( Finteraction, as.numeric(ag[[3]][4]))
      
      # derive p values using the pf function, to avoid zero p values of which
      # you can not take a log
      pmarker      <- c( pmarker,      pf(1/as.numeric(ag[3,3]),ag[3,2],ag[3,1]) )
      pinteraction <- c( pinteraction, pf(1/as.numeric(ag[4,3]),ag[4,2],ag[4,1]) )    
    }
  }

  # collect results
  perprobeset <- data.frame( markers=rownames(genotypes), 
                             Fmarker=Fmarker, Finteraction=Finteraction, 
                             pmarker=pmarker, pinteraction=pinteraction )

  if (!is.null(filename))    
  {
    write.table( perprobeset, file=filename, sep="," )
  }
                             
  return( perprobeset ) #return results
}
