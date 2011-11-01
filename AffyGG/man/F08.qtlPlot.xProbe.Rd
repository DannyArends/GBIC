\name{qtlPlot.xProbe}
\alias{qtlPlot.xProbe}

\title{ QTL plots of Affymetrix probe level data }

\description{
  QTL plots of Affymetrix probe level data.
}

\usage{

  qtlPlot.xProbe ( probesetName, markersPos, probeQtlProfiles,
                   qtlThres, chrOffsets, filename=NULL ) 
}

\arguments{
  \item{probesetName}{ name of the probe set. }
  \item{markersPos}{ marker positions in basepairs (bp). }
  \item{probeQtlProfiles}{ data.frame with markers, probe numbers and
                           \eqn{-log_{10}(p)}{-log10(p)} values of marker effect. }
  \item{qtlThres}{ threshold value for marker effect. }
  \item{chrOffsets}{ starting positions per chromosome in basepairs. }  
  \item{filename}{ (optional) path of the file where the data will be loaded from or save to.\cr
   In this function, png output file where results plot is written. Default is NULL. } 
}

\details{
  \code{\link{qtlPlot.xProbe}} plots a genome wide QTL profile for each
  individual probe.
}

\value{
  It generates a plot (or PNG image file) with the QTL profile of all probes within
  one probe set.
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
    \code{\link{qtlPlot.xProbeSet}}, \code{\link{qtlMap.xProbe}}
}
  
\examples{
# We load the QTL profile of our probes
data( probeQtlProfiles )

# We calculate the -log(p) values on the pmarker column
qtlProfiles <- cbind( probeQtlProfiles[,1:3], minlog10pmarker=-log10(probeQtlProfiles[,4]))

# We load the positions of our markers
data( markersPos )

# We indicate the starting positions of the chromosomes in bp
data( chrOffsets ) 

# We create a plot with the QTL profiles of the probes of one probe set 
qtlPlot.xProbe( "96243_f_at", markersPos, qtlProfiles, 
                3.72, chrOffsets )
}