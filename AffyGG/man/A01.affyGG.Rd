\name{affyGG}
\alias{affyGG}
\title{ Protocol function to analyze Affymetrix genetical genomics data }
\description{
  This function provides a generic protocol to analyze genetical genomics data
  obtained from experiments done with Affymetrix microarrays. The protocol
  includes pre-processing of the Affymetrix data, QTL analysis, significance threshold calculation and
  visualization of probe signals and QTL profiles. Finally, it includes a
  procedure to backwards eliminate deviating probes from probe sets.                  
}
\usage{
affyGG( genotypes, batch, nProbes, nIndiv, n.simulations, 
        directory, celfiles, probelevel, probesetName, 
        overwrite, n.Plots, markersPos, probesPos, chrOffsets)
}

\arguments{
  \item{genotypes}{ matrix of genotypes for each marker (rownames) and individual 
                    (columnnames), as numeric values 1, 2 or NA when missing.\cr
                    See \code{\link{genotypes}} example data.}
  \item{batch}{ vector with batch numbers in which each individual (e.g. BXD)
                was processed. }
  \item{nProbes}{ (optional) number of probes per probeset. Default is 16. }
  \item{nIndiv}{ (optional) number of individuals. Default is 30. }
  \item{n.simulations}{ (optional) number of simulations. Default is 100 times. }
  \item{directory}{ path of the directory where files are located.\cr
                    In this function, input directory where .CEL files are located. }
  \item{celfiles}{ name of the .CEL files that should be processed. }
  \item{probelevel}{ (optional) boolean value to indicate whether the analysis 
                     should be done at probe level (\code{T}) or at probe set 
                     level (\code{F}). Default is \code{T}. }
  \item{probesetName}{ (optional) name of the probe set. Default is NULL. }
  \item{overwrite}{ (optional) boolean value to indicate whether the files 
                    with previous results are kept (\code{F}) or overwriten (\code{T}).
                    Default is \code{F}.  }
  \item{n.Plots}{ (optional) number of plots to be created when a probe set name (i.e. \code{probesetName})
                  is not indicated. Default is 5. }
  \item{markersPos}{ marker positions in basepairs (bp). }                  
  \item{probesPos}{ interrogation positions of the probes on the mRNA. }  
  \item{chrOffsets}{ starting positions per chromosome in basepairs. }  
}

\details{
First, affyGG perform pre-processing of the Affymetrix data in the form of RMA
background correction and quantile normalization (Irizarry et al. 2003). Next,
QTL analysis is performed by a statistical multi-probe model introduced by
Alberts et al. (2005). This QTL analysis can be performed on probe set level,
obtaining significance level for QTL effect and QTL \eqn{\times}{x} probe interaction effect.
Also, QTL analysis can be performed on individual probes separately, obtaining
QTL profiles per probe. If microarrays were processed in
batches, a factor batch can be included in the statistical model to correct
for batch effects in the data. By simulation significance thresholds for QTL
effects are calculated. The probe signals of individual probes within one
probe set are plotted, including the relative probe positions on the
mRNA. Finally, QTL profiles for the analyses on probe set level and on
probe level are plotted along the genome. 
}
\value{
  \item{probesignals}{Create a .csv file where on each line the probe signal values for all .CEL files supplied are
  given.}
  \item{qtlThres}{Returns a .csv file where on each 
  line (as many lines as \code{n.simulations}) the highest F and p values for marker 
  and interaction effect for a genome scan is given. The 95th percentile is
  chosen as threshold.}
  \item{qtlMapProbeset}{ Returns a file containing QTL mapping results on
  probe set level: for each marker, 
  the F-value and p-value of both the marker and the marker x probe effect (interaction) are given. }
  \item{qtlMapProbe}{  Returns a file containing QTL mapping results on probe level: for each marker, 
  the F-value and p-value of the marker effect is given.}
  \item{probePlot}{ Generates a .png file with the probe signals of individual
  probes for one probe set
  with a different color per allele (if supplied).}
  \item{qtlPlotProbeset}{ Generates a .png file with the QTL profile on probe set level
  and a profile of the interaction effect.}
  \item{qtlPlotProbe}{ Generates a .png file with the QTL profile of all probes within
  one probe set.}
  Example:  
  \preformatted{
Package: affyGG
Computational protocols for genetical genomics analyses with 
Affymetrix arrays.
                                                                               
STEP 1: Checking arguments... (T: 0.01, P: 0.01 secs elapsed)             
STEP 2: Preprocessing Affymetrix .CEL files... (T: 86.62, 
        P: 86.61 secs elapsed)
STEP 3.1: Significance thresholds calculation... (T: 221.62, 
        P: 135 secs elapsed)
STEP 3.2: Choosing the 95 percentile value... (T: 221.62, 
        P: 0 secs elapsed) 
STEP 4.1: QTL analysis... (T: 235.16, P: 13.54 secs elapsed)                 
STEP 5: Plotting probe signals... (T: 235.44, P: 0.28 secs elapsed)          
STEP 6: Plotting QTL maps... (T: 235.59, P: 0.15 secs elapsed)               
  
  }
}

\references{ 
  Alberts R, Terpstra P, Bystrykh LV, de Haan G, Jansen RC. A statistical
  multiprobe model for analyzing cis and trans genes in genetical genomics
  experiments with short-oligonucleotide arrays. Genetics. 2005
  Nov;171(3):1437-9. \cr
  
  Irizarry RA, Hobbs B, Collin F, Beazer-Barclay YD,
  Antonellis KJ, Scherf U, Speed TP. Exploration, normalization, and
  summaries of high density oligonucleotide array probe level data. 
  Biostatistics. 2003 Apr;4(2):249-64.
}

\author{ 
Rudi Alberts <r.alberts@rug.nl>, \cr
Gonzalo Vera <g.vera.rodriguez@rug.nl>, \cr
Ritsert Jansen <r.c.jansen@rug.nl>
}

\seealso{ 
  \code{\link{rma.preprocessing}}, \code{\link{qtlThresholds.sma}},\cr 
  \code{\link{qtlMap.xProbe}}, \code{\link{qtlMap.xProbeSet}}, \cr
  \code{\link{probePlot}},\cr
  \code{\link{qtlPlot.xProbe}}, \code{\link{qtlPlot.xProbeSet}}
}
\examples{
\dontrun{
# We load and set our arguments variables
data(genotypes)
data(batch)
data(markersPos)
data(chrOffsets)

directory <- "c:/data/my_project"

celfiles <- c( "bxd1a.cel", "bxd2a.cel", "bxd5a.cel", "bxd6a.cel", 
               "bxd8a.cel", "bxd9a.cel", "bxd11a.cel", "bxd12a.cel", 
               "bxd14a.cel", "bxd15a.cel", "bxd16a.cel", "bxd18a.cel", 
               "bxd19a.cel", "bxd21a.cel", "bxd22a.cel", "bxd24a.cel", 
               "bxd25a.cel", "bxd27a.cel", "bxd28a.cel", "bxd29a.cel", 
               "bxd30a.cel", "bxd31a.cel", "bxd32a.cel", "bxd33a.cel", 
               "bxd34a.cel", "bxd36a.cel", "bxd38a.cel", "bxd39a.cel", 
               "bxd40a.cel", "bxd42a.cel" )
               
probesPos <- c( 1893, 1894, 1897, 1904, 1906, 1911, 1912, 1913, 
                1916, 1925, 1929, 1933, 1934, 1944, 1948, 1955)


# If we don't provide a probe set name, it will process all the probe sets,
# which it takes several hours, although it is only needed to process once.
probesetName <- "96243_f_at"                          
                               
# For testing we decrease the number of markers and simulations, better 
# than the number of individuals, to decrease processing time. 
affyGG( genotypes=genotypes[1:100,], batch=batch, directory=directory,
        celfiles=celfiles, probesetName=probesetName, n.simulations = 20,
        markersPos=markersPos[1:100], probesPos=probesPos, 
        chrOffsets=chrOffsets)

}
}
