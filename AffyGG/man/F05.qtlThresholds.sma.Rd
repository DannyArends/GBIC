\name{qtlThresholds.sma}
\alias{qtlThresholds.sma}

\title{ Determine QTL threshold via simulation of expression data }

\description{
  Determine QTL threshold using simulated expression data and the real 
  genotype data.
}

\usage{
  qtlThresholds.sma( genotypes, batch, 
                     nProbes=16, nIndiv=30, 
                     n.simulations=10000, filename=NULL )
}

\arguments{

  \item{genotypes}{ matrix of genotypes for each marker (rownames) and individual 
  (columnnames), as numeric values 1, 2 or NA when missing.\cr
  See \code{\link{genotypes}} example data.}

  \item{batch}{ vector with the batch number in which each individual (e.g. BXD)
  was processed.}

  \item{nProbes}{ (optional) number of probes per probeset. Default is 16.}
  
  \item{nIndiv}{ (optional) number of individuals. Default is 30.}
  
  \item{n.simulations}{ (optional) number of simulations. Default is 10000 times. }
  
  \item{filename}{ (optional) path of the file where the data will be loaded from or save to.\cr
  In this function, output file where results are stored (*.csv = comma
  separated values). }       
}

\details{
  Random probe level data is generated (mean 0 and sd 1) and
  with the real markers, QTL analyses are performed 10000 times. The
  highest QTL value per genome scan is stored. Finally the 95\% value can
  be used as threshold.
}

\value{
  Returns a \code{data.frame} (and optionally creates a .csv file) where on each 
  line (as many lines as \code{n.simulations}) the highest F value for marker 
  and interaction effect for a genome scan is given. \cr
  Example of 5 simulations: \cr
\tabular{lllll}{
    \tab maxFmarker \tab maxFinteraction \tab maxpmarker \tab maxpinteraction \cr
1   \tab 8.607366        \tab 2.376508  \tab 0.9999094       \tab 0.9957732   \cr
2   \tab 6.599920        \tab 2.868144  \tab 0.9983618       \tab 0.9989851   \cr
3   \tab 6.070313        \tab 2.439586  \tab 0.9987627       \tab 0.9991177 \cr
4  \tab 10.392663        \tab 2.649089  \tab 0.9985993       \tab 0.9994769 \cr
5   \tab 5.831100        \tab 2.493569  \tab 0.9964620       \tab 0.9999437 \cr
}
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
   \code{\link{qtlMap.xProbeSet}},  \code{\link{qtlMap.xProbe}}
}
\examples{
# We load the genotypes and batch numbers from the included example data
data(genotypes)
data(batch)

# We execute the simulation. Around 20 minutes per simulation on Pentium 4
qtlTh <- qtlThresholds.sma( genotypes, batch, 16, 30, 5 )

head( qtlTh )
}