\name{qtlPlot.xProbeSet}
\alias{qtlPlot.xProbeSet}
\title{ QTL plots of Affymetrix probe set level data }

\description{
  QTL plots of Affymetrix probe set level data. 
}

\usage{                      
   qtlPlot.xProbeSet( probesetName, markersPos, 
                      probesetQtlProfile, interactionProfile,
                      qtlThres, interactionThres, 
                      chrOffsets, filename=NULL ) 
}

\arguments{
  \item{probesetName}{ name of the probe set. }
  \item{markersPos}{ marker positions in basepairs (bp). }
  \item{probesetQtlProfile}{ data.frame with markers and
        \eqn{-log_{10}(p)}{-log10(p)} values of marker effect. }
  \item{interactionProfile}{ data.frame with markers, and
         \eqn{-log_{10}(p)}{-log10(p)} values of interaction effect between marker and probe. }
  \item{qtlThres}{ threshold value for marker effect. }
  \item{interactionThres}{ threshold value for interaction effect. }
  \item{chrOffsets}{ starting positions per chromosome in basepairs. }  
  \item{filename}{ (optional) path of the file where the data will be loaded from or save to.\cr
   In this function, png output file where results plot is written. Default is NULL. } 
}

\details{                                                  
  \code{\link{qtlPlot.xProbeSet}} plots a genome wide QTL profile, with marker
  effect in blue and marker by probe interaction in green.   
}

\value{
  It generates a plot (or PNG image file) with the QTL profile on probe set level
  and a profile of the interaction effect.
}

\references{
  Alberts R, Terpstra P, Bystrykh LV, de Haan G, Jansen RC. A statistical 
  multiprobe model for analyzing cis and trans genes in genetical genomics 
  experiments with short-oligonucleotide arrays. Genetics. 2005 Nov;171(3):1437-9.
}

\author{ 
Rudi Alberts <r.alberts@rug.nl>, 
Ritsert Jansen <r.c.jansen@rug.nl>
}

\seealso{
    \code{\link{qtlPlot.xProbe}}, \code{\link{qtlMap.xProbeSet}}
}
\examples{
# We retrieve the output of gg.qtlMap.xProbeSet
data( probesetQtlProfiles )

# We calculate the -log10 of the marker and interaction p-values
interactionProfile <- -log10( probesetQtlProfiles$pinteraction )
qtlProfile <- -log10( probesetQtlProfiles$pmarker )

# We load the positions of our markers
data( markersPos )

# We indicate the starting positions of the chromosomes in bp
data( chrOffsets )
                  
# We plot the QTL profile on probe set level and a profile of the 
#   interaction effect.
qtlPlot.xProbeSet( "96243_f_at", markersPos, qtlProfile, 
                   interactionProfile, 3.65, 3.75, chrOffsets )
}