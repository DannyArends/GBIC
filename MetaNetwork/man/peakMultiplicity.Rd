\name{findPeakMultiplicity}
\title{Cluster multiple mass peaks based on correlation and mass}
\alias{findPeakMultiplicity}
\description{ Predict metabolites that may be represented in multiple mass peaks 
based on correlation and m/z value. If two mass peaks are highly correlated 
(r>0.95) and their mass difference 
is 1 or 2 (as isotopes) or mass ratio is 2, 3, or 1/2, 1/3 (as different charged), 
they are predicted as multiple peaks of a same metabolite.}
\usage{
findPeakMultiplicity(corrZeroOrder, peaks, corrThres=0.95, filename=NULL)
}
\arguments{
  \item{corrZeroOrder}{ a zero order peak correlation matrix that has been calculated before.
  See \code{\link{corrZeroOrder}} example data.}
  \item{peaks}{ matrix of mass/charge peaks (column1) for each trait (rownames). 
  See \code{\link{peaks2}} example peaks data for unidentified example traits data \code{\link{traits2}}.}
  \item{corrThres}{ (optional) numeric threshold for significant correlated peaks
  that should be tested for multiple peaks relationships. Default is 0.95. }
  \item{filename}{ (optional) path of the csv file where the multiple peak summary 
  is to be stored. Default NULL. }
}
\value{
  Returns a data frame with Peak Multiplicity summary containing the following headers:
  \item{cluster}{ cluster id of related multiple peaks.}   
  \item{peak1}{ mass over charge peak.}    
  \item{mz1}{ mass over charge ratio of \code{peak1}.}
  \item{peak2}{ mass over charge peak that is related to \code{peak1}.}
  \item{mz2}{ mass over charge ratio of \code{peak2}. }    
  \item{corrCoef}{ correlation coefficient between \code{peak1} and \code{peak2}. }
  \item{massDiff}{ peak difference: \code{mz1} - \code{mz2}.}
  \item{massRatio}{ peak ratio: \code{mz1}/\code{mz2}.}
  \item{relationship}{ putative \code{isotope} or \code{diffCharged} relationship.}
  See \code{\link{peakMultiplicity}} example data set.
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
\note{
The names of traits should be consistent over \code{corrZeroOrder} and \code{peaks}.
}
\seealso{ 
Use \code{\link{traits2}} and \code{\link{peaks2}} as example unidentified peak 
data set or use \code{\link{loadData}} to load your own data.
Use \code{\link{qtlCorrZeroOrder}} to calculate \code{corrZeroOrder} zero order correlation.
Use \code{\link{MetaNetwork}} for automated application of this function as part a genetic analysis protocol on metabolites.
}
\examples{
## load the example data provided with this package                         
data(traits2)   
data(peaks2)    
data(genotypes)     
data(markers)                                         
                                             
##OR: load your own data                     
#traits2      <- loadData("traits2.csv")
#peaks2       <- loadData("peaks2.csv")
#genotypes    <- loadData("genotypes.csv")
#markers      <- loadData("markers.csv")
                                             
##calculate the two part qtl
qtlProfiles   <- qtlMapTwoPart(genotypes=genotypes, traits=traits2, spike=4)
  
##set the qtl threshold
qtlThres      <- 3.79

##OR: estimate the threshold yourself
#qtlThres     <- qtlThreshold(genotypes, traits, spike=4)
  
##calculate zero order correlation
corrZeroOrder <- qtlCorrZeroOrder(markers, qtlProfiles, qtlThres)

##find peak multiplicity
peakMultiplicity <- findPeakMultiplicity(corrZeroOrder, peaks2)

##show summary
peakMultiplicity
}