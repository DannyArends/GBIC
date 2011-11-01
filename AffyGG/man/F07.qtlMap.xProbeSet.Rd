\name{qtlMap.xProbeSet}
\alias{qtlMap.xProbeSet}
\title{QTL mapping of Affymetrix probe set level data}

\description{
  Probe set level QTL analysis of Affymetrix data.
}

\usage{
  qtlMap.xProbeSet( genotypes, traits, batch, filename )
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
  \code{\link{qtlMap.xProbeSet}} fits a mixed effects model per probe set
  with fixed effects batch (optional), probe, batch x probe (optional), marker, marker x probe and
  random effect individual.\cr
}

\value{
  Returns a \code{data.frame} containing QTL mapping results: for each marker, 
  the F-value and p-value of both the marker and the marker x probe effect (interaction) are given. \cr
  Example:
  \preformatted{
> qtlMap <- qtlMap.xProbeSet( genotypes, traits, batch )
> qtlMap[1:5,]}\cr

  \tabular{lrrrrr}{  
      \tab markers    \tab Fmarker \tab Finteraction   \tab pmarker \tab pinteraction \cr
1 \tab D01Msw003 \tab 0.00778279   \tab 0.06869318 \tab 0.9303777    \tab 0.9999997 \cr
2 \tab D01Msw011 \tab 1.08072774   \tab 1.15390195 \tab 0.3081063    \tab 0.3061359 \cr
3 \tab D1Mit294 \tab 1.08072774   \tab 1.15390195 \tab 0.3081063    \tab 0.3061359 \cr
4 \tab D1Mit67 \tab 0.12549663   \tab 0.11241176 \tab 0.7260048    \tab 0.9999896 \cr
5 \tab D1Mit430 \tab 0.04364617   \tab 0.12714113 \tab 0.8361415    \tab 0.9999763 \cr
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
Rudi Alberts <r.alberts@rug.nl>, \cr
Gonzalo Vera <g.vera.rodriguez@rug.nl>, \cr
Ritsert Jansen <r.c.jansen@rug.nl>
}

\seealso{
  \code{\link{qtlMap.xProbe}}

}
\examples{
# We load the genotypes, traits and batch numbers from the example data
data(genotypes)
data(traits)
data(batch)

# We perform QTL mapping of one probe set and show the first results
qtlMap <- qtlMap.xProbeSet( genotypes, traits, batch )
head( qtlMap, 20 )
}