\name{peakMultiplicity}
\docType{data}
\alias{peakMultiplicity}
\title{Example Peak Multiplicity report}
\description{
Example of the Multiple Peaks report that relates metabolites represented 
by multiple mass peaks (because of isotopic or differential charges) for
\code{\link{peaks2}} example data set using information from \code{\link{traits2}}.
This report can be produced 
by functions \code{\link{findPeakMultiplicity}} or \code{\link{MetaNetwork}}.
\preformatted{
> data(peakMultiplicity)
> peakMultiplicity[1:5,]}
\tabular{llrrrrrrr}{
cluster \tab peak1       \tab mz1 \tab peak2        \tab mz2 \tab corrCoef  \tab massDiff \tab massRatio \tab relationship \cr
     1  \tab LCavg.1539  \tab 378 \tab LCavg.1538   \tab 377 \tab 0.9998261 \tab 1        \tab 1.0       \tab isotope      \cr
     2  \tab LCavg.1570  \tab 424 \tab LCavg.1566   \tab 423 \tab 0.9989359 \tab 1        \tab 1.0       \tab isotope      \cr
     3  \tab LCavg.1610  \tab 436 \tab LCavg.1596   \tab 873 \tab 0.9912291 \tab -437     \tab 0.5       \tab diffCharged  \cr
     4  \tab LCavg.1611  \tab 437 \tab LCavg.1610   \tab 436 \tab 0.9612162 \tab 1        \tab 1.0       \tab isotope      \cr
     5  \tab LCavg.1612  \tab 438 \tab LCavg.1596   \tab 873 \tab 0.9513867 \tab -435     \tab 0.5       \tab diffCharged  \cr
}

\code{peakMultiplicity} has the following headers:
\itemize{
  \item{cluster}{ cluster id of related multiple peaks.}   
  \item{peak1}{ mass over charge peak.}    
  \item{mz1}{ mass over charge ratio of \code{peak1}.}
  \item{peak2}{ mass over charge peak that is related to \code{peak1}.}
  \item{mz2}{ mass over charge ratio of \code{peak2}. }    
  \item{corrCoef}{ correlation coefficient between \code{peak1} and \code{peak2}. }
  \item{massDiff}{ peak difference: \code{mz1} - \code{mz2}.}
  \item{massRatio}{ peak ratio: \code{mz1}/\code{mz2}.}
  \item{relationship}{ putative \code{isotope} or \code{charge} relationship.}
}
}
\usage{data(peakMultiplicity)}
\format{ 
  A data frame.
}
\source{ 
  Keurentjes JJB, FU J, de vos CHR, Lommen A, Hall RD, Bino RJ, van der Plas LHW, Jansen RC, Vreugdenhil D
  and Koornneef M. The genetics of plant metabolism. Nature Genetics (2006) 7: 842-849. 
}
\references{  
  Fu J, Swertz MA, Keurentjes JJB, Jansen RC. MetaNetwork: a computational tool for the 
  genetic study of metabolism. Nature Protocols (2007). 
  
  \url{http://gbic.biol.rug.nl/supplementary/2007/MetaNetwork}
}
\author{
Jingyuan Fu <j.fu@rug.nl>, 
Morris Swertz <m.a.swertz@rug.nl>, 
Ritsert Jansen <r.c.jansen@rug.nl>
}
\seealso{ 
  Use \code{\link{findPeakMultiplicity}} to generate \code{peakMultiplicity}.
  Use \code{\link{MetaNetwork}} to produce this dataset as part a genetic 
  analysis protocol on metabolites.
}
\examples{
## load the example data provided with this package
data(peakMultiplicity)

##OR: generate your own qtlSummary using MetaNetwork
#genotypes  <- loadData("genotypes.csv") 
#traits2    <- loadData("traits2.csv") 
#markers    <- loadData("markers.csv")  
#peaks2     <- loadData("peaks2.csv")
#
#MetaNetwork (markers=markers, genotypes=genotypes, traits=traits2, peaks=peaks2,
#             spike=4) 
             
##view the summary
peakMultiplicity              
}
\keyword{datasets}