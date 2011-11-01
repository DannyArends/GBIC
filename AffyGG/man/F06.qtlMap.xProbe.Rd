\name{qtlMap.xProbe}
\alias{qtlMap.xProbe}
\title{ QTL mapping of Affymetrix probe level data }

\description{
  Probe level QTL analysis of Affymetrix data.
}

\usage{
  qtlMap.xProbe( genotypes, traits, batch, filename )
}

\arguments{

  \item{genotypes}{ matrix of genotypes for each marker (rownames) and individual 
  (columnnames), as numeric values 1, 2 or NA when missing.\cr
  See \code{\link{genotypes}} example data.}
  
  \item{traits}{ matrix of phenotypes for each trait (rownames) and individual 
  (columnnames), as numeric or NA when missing.\cr
  See \code{\link{traits}} example data. }

  \item{batch}{ (optional) vector with batch numbers in which each individual (e.g. BXD)
  was processed. Default is NULL.}

  \item{filename}{ (optional) path of the file where the data will be loaded from or save to.\cr
  In this function, output file where results are stored (*.csv = comma
  separated values). Default is NULL. } 
}

\details{
  \code{\link{qtlMap.xProbe}} fits a fixed effects model per probe with
  factor batch (optional) and factor marker.
}

\value{
  Returns a \code{data.frame} containing QTL mapping results: for each marker, 
  the F-value and p-value of the marker effect is given.\cr
  Example:
\preformatted{
> qtlMap <- qtlMap.xProbe( genotypes, traits, batch )
> qtlMap[1:5,]}\cr
\tabular{lllll}{
       \tab marker \tab probenr \tab Fmarkerperprobe \tab pmarkerperprobe \cr
1 \tab D01Msw003       \tab 1    \tab 0.5021788295      \tab 0.48485016 \cr
2 \tab D01Msw003       \tab 2    \tab 0.1797457271      \tab 0.67507905 \cr
3 \tab D01Msw003       \tab 3    \tab 3.9582922214      \tab 0.05725642 \cr
4 \tab D01Msw003       \tab 4    \tab 0.0318835120      \tab 0.85966763 \cr
5 \tab D01Msw003       \tab 5    \tab 0.0002685086      \tab 0.98705136 \cr  
}  
}

\references{
  Alberts R, Terpstra P, Bystrykh LV, de Haan G, Jansen RC. A statistical
  multiprobe model for analyzing cis and trans genes in genetical genomics
  experiments with short-oligonucleotide arrays. Genetics. 2005 Nov;171(3):1437-9.

}

\note{
  Argument \code{traits} expects a data.frame with probe signal values 
  of \eqn{one} probe set for each individual (columnnames) and each probe (rownames).
}

\author{
Rudi Alberts <r.alberts@rug.nl>,   \cr
Gonzalo Vera <g.vera.rodriguez@rug.nl>,\cr
Ritsert Jansen <r.c.jansen@rug.nl>
}

\seealso{
  \code{\link{qtlMap.xProbeSet}}

}
\examples{
# We load the genotypes, traits and batch numbers from the example data
data(genotypes)
data(traits)
data(batch)

# We perform QTL mapping of one probe set and show the first results
qtlMap <- qtlMap.xProbe( genotypes, traits, batch )
head( qtlMap, 20 )
}
