# Package: GeneticalGenomics v1.0.0
# Originally part of Rudi Alberts' PhD Thesis
# Function name: gg.probeElimination
# Author: Rudi Alberts <r.alberts@rug.nl>
# Maintainer: Gonzalo Vera <g.vera.rodriguez@rug.nl>
# Version: 1.0.0
# Date: 1 Feb. 2006

probeElimination <- function( probesetName, markerName, genotypes, traits, 
                                 batch, filename=NULL, appendRecord=F ) 
{
  if( is.null(filename) ){
      filename <- paste( "drops.", auxCleanFilename( probesetName ), 
                         ".csv", sep="" )
  }       
  genotype <- genotypes[ markerName, ]

  nProbes <- nrow( traits )           # nr. of probes        
  onevec  <- c( as.matrix( traits ) ) # put data in one vector
  
  anovain <- data.frame( batch= rep( batch, each=nProbes ), 
                         probe= 1:nProbes,
                         mouse= rep( 1:ncol( traits ), each=nProbes ), 
                         value= onevec,
                         marker = rep( genotype, each=nProbes) )
  
  interactionEffect <- dropprobes( -1, anovain )[1,2]
  
  cat( probesetName, ",,,, startint: ,", interactionEffect, "\n", 
       file=filename, append=appendRecord )
  
  results <- data.frame( 
                  dropped= 0, 
                  pmarker= 0,
                  pinteraction= interactionEffect,
                  pmarkerdropped= 0, 
                  pintdropped= 0, 
                  alarm= F,
                  allele1bigger= F, 
                  allele1biggerdropped= F )

  ## probe backwards elimination procedure: 
  ##  drop probes one by one, 
  ##  see for which probe the interaction effect is reduced most, 
  ##  leave this probe out permanently (it is added to 'dropped'). 
  ##  repeat for the remaining probes. 
  ##  store qtl and interaction effect after each drop.
                             
  dropped <- NULL
  # drop until 2 probes left and store all results
  while( length(dropped) < (nProbes-2) ) 
  { 
    newint <- NULL
    # drop each probe one by one that is not yet permanently dropped
    for( k in 1:nProbes ) 
    { 
      if( !k %in% dropped ) 
      { 
        doit    <- dropprobes( c(k,dropped), anovain )
        newint  <- rbind( newint, 
                          data.frame( probe=k, marker=doit[1,1], int=doit[1,2],
                                      allele1bigger=doit[1,3]))
      }
    }
    # determine for which probe the interaction effect became the lowest.
    newint  <- newint[order(newint$int,decreasing=T),] 
    # store this probe in 'dropped'
    dropped <- c(dropped,newint[1,1]) 
    int     <- newint[1,3]

    ## what happened to the dropped probes - here drop the non-dropped 
    ##  probes to determine this 
        
    qtldropped           <- intdropped  <- -1
    allele1biggerdropped <- F
    restant              <- 1:nProbes
    restant              <- restant[ !restant %in% dropped ]
    
    # At least we need 2 probes in order to calculate the interaction effect
    if( length(dropped) > 1 ) 
    {
      doit                 <- dropprobes( restant, anovain )
      qtldropped           <- doit[1,1]
      intdropped           <- doit[1,2]
      allele1biggerdropped <- doit[1,3]
    } 

    # store results in a file per probe set
    cat( "dropped: ,",    newint[1,1], 
         ",pmarker: ,",    newint[1,2],
         ",pinteraction: ,",       newint[1,3],
         ",pmarkerdropped: ,", qtldropped, 
         ",pintdropped: ,", intdropped, 
         ",alarm: ,",     newint[1,2] > qtldropped,
         ",allele1bigger: , ",  newint[1,4], 
         ",allele1biggerdropped: ,", allele1biggerdropped,
         "\n", file=filename, append=T )
      
      results<-rbind(results,data.frame( 
                  dropped= newint[1,1], 
                  pmarker= newint[1,2],
                  pinteraction= newint[1,3],
                  pmarkerdropped= qtldropped, 
                  pintdropped= intdropped, 
                  alarm= (newint[1,2] > qtldropped),
                  allele1bigger= newint[1,4], 
                  allele1biggerdropped= allele1biggerdropped ))
  }
  
}

# todrop like c(1,5,7,8,9) (indexes)
# todrop -> which probes we want to exclude
# anovain -> data.frame 
dropprobes <- function( todrop, anovain ) 
{
  library(nlme)
  # needed to make parameter estimates add up to zero
  options( contrasts=c("contr.sum","contr.poly") )
  
  temp          <- anovain[ ! anovain$probe %in% todrop,]
                   # fit linear mixed-effects model
  g             <- lme( value~factor(batch) + factor(probe) + 
                        factor(batch):factor(probe) + 
                        factor(marker) + factor(marker):factor(probe),
                        data=temp, random=~1|factor(mouse) )  
  ag            <- anova(g)  # compute anova table
  pmarker       <- pf( 1/ag[[3]][4], ag[[2]][4], ag[[1]][4] )
  pinteraction  <- pf( 1/ag[[3]][6], ag[[2]][6], ag[[1]][6] )
  
  # chech whether B6 has higher (1) or lower (0) expression values than D2
  allele1bigger      <-  
    mean( temp[temp$marker==1,]$value ) > mean(temp[temp$marker==2,]$value)

  return( data.frame(pmarker, pinteraction, allele1bigger) )
}

