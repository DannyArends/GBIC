\name{probePlot}
\alias{probePlot}
\title{ Plots of Affymetrix probe level data }

\description{ Plots of Affymetrix probe level data. }

\usage{
probePlot( traits, probesetName, probesPos, alleleColors, filename ) 
}

\arguments{
  \item{traits}{ matrix of phenotypes for each trait (rownames) and individual 
  (columnnames), as numeric or NA when missing.\cr
  See \code{\link{traits}} example data. }
  
  \item{probesetName}{ name of the probe set. }

  \item{probesPos}{ interrogation positions of the probes on the mRNA. }  

  \item{alleleColors}{ (optional) vector with R color codes indicating the allele that 
  each individual is carrying for one marker. Default is 1 (black).}

  \item{filename}{ (optional) path of the file where the data will be loaded from or save to.\cr
  In this function, png output file where results plot is written. Default is NULL. }
}

\details{
  \code{\link{probePlot}} plots Affymetrix probe signals per probe set for
  multiple chips. Any type of signals (raw, pre-processed, PM, MM) can
  be input.
}

\value{
  It generates a plot (or PNG image file) with the probe signals for one probe set
  with a different color per allele (if supplied).                                                    
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
  \code{\link{rma.preprocessing}}
}
\examples{
# We indicate the bp positions of the probes on the mRNA
pos<-c(1893, 1894, 1897, 1904, 1906, 1911, 1912, 1913, 1916, 1925, 1929, 1933, 1934, 1944, 1948, 1955)

# We retrieve the different alleles for one marker 
# (these values: 1 or 2, can be used like color codes)
data( genotypes )
alleleColors <- as.numeric( genotypes["D1Mit145",] )

# We load the traits of the probeset
data( traits ) 

# We plot the probe signals
probePlot( traits, "96243_f_at", pos, alleleColors )

}
