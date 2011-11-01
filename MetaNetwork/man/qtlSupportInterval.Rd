\name{qtlSupportInterval}
\alias{qtlSupportInterval}
\title{ Estimate QTL support interval }
\description{
  For one QTL profile, significant QTLs are selected based on \code{qtlThres}. 
  The regions within the \code{interval.dropoff} of these QTL peaks
  are defined as a support interval. The most left and the most right marker within
  each support interval are returned as matrix.
}
\usage{
qtlSupportInterval(markers, oneQtlProfile, qtlThres, interval.dropoff = 1.5)
}
\arguments{
  \item{markers}{ matrix of markers (rownames) and their chromosome numbers 
  (column 1) and centi-Morgan positions (cM, column 2), ordered by position.
  See \code{\link{markers}} example data.}
  \item{oneQtlProfile}{ one row from the \code{\link{qtlProfiles}} matrix of QTL mapping of traits (rownames) 
  to markers (columnnames), as \eqn{-log_{10}(p)} values.
  See \code{\link{qtlProfiles}} example data.}
  \item{qtlThres}{ numeric \eqn{-log_{10}(p)} threshold value for significant QTLs. }
  \item{interval.dropoff}{ numeric drop-off \eqn{-log_{10}(p)} value from the QTL peak that 
  defines left and right border of support interval (optional). Default is 1.5. }
}
\value{
  Returns a matrix of the markers number that define the left and right borders of each
  support interval for \code{oneQtlProfile}.
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
The markers should be ordered sequentially. The names of markers should be 
consistent over \code{markers} and \code{qtlProfiles}.
}
\seealso{ 
Use \code{\link{markers}} as example data set or use \code{\link{loadData}} 
to load your own data.
Use \code{\link{qtlMapTwoPart}} for the calculation of \code{qtlProfiles}. 
Use \code{\link{qtlThreshold}} for the estimation of \code{qtlThres} QTL 
significance threshold.
Use \code{\link{qtlSummary}} for automated application of this function to 
produce a support interval summary for a set \code{qtlProfiles}.
Use \code{\link{MetaNetwork}} for automated appliction of this function as 
part a genetic analysis protocol on metabolites.
}
\examples{
## load the example data provided with this package                          
data(markers)   
data(genotypes)
data(traits)                                                       
                                             
##OR: load your own data                     
#markers    <- loadData("markers.csv")
#genotypes  <- loadData("genotypes.csv")
#traits     <- loadData("traits.csv")  
                                             
##calculate the two part qtl for only the first trait
qtlProfiles <- qtlMapTwoPart(genotypes=genotypes, traits=traits[1,], spike=4)

##set the qtl threshold
qtlThres    <- 3.79

##OR: estimate the threshold yourself
#qtlThres   <- qtlThreshold(genotypes, traits, spike=4)

##calculate qtl confidence interval for the first qtl profile
qtlSuppInt  <- qtlSupportInterval(markers, qtlProfiles[1,], qtlThres=qtlThres, 
                                  interval.dropoff = 1.5)    
                                  
##show the interval
qtlSuppInt                                                             
}