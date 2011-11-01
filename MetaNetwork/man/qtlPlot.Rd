\name{qtlPlot}
\alias{qtlPlot}
\title{ Plot QTL profiles  }
\description{
  Plot QTL profiles. 
}
\usage{
qtlPlot(markers, qtlProfiles, qtlThres, addTitle=NULL, addMarkerLabels=FALSE, addLegend=TRUE, color=NULL)
}
\arguments{
  \item{markers}{ matrix of markers (rownames) and their chromosome numbers (column 1) 
  and centi-Morgan (cM) positions, ordered by position.
  See \code{\link{markers}} example data.}
  \item{qtlProfiles}{ matrix of QTL mapping of traits (rownames) 
  to markers (columnnames), as \eqn{-log_{10}(p)} values.
  See \code{\link{qtlProfiles}} example data.}
  \item{qtlThres}{ \eqn{-log_{10}(p)} threshold value for significant QTLs. }
  \item{addTitle}{ (optional) show a title for the plot. Default is NULL. }
  \item{addMarkerLabels}{ (optional) show markernames as labels. Default is FALSE. }
  \item{addLegend}{ (optional) show a legend. Default is TRUE. }
  \item{color}{ (optional) a vector of colors for the plot such as \code{c(1,2,3)} 
  (when shorter than length of QTLprofiles the colors will be applied repeatedly). Default NULL. 
}
}
\value{
  Draws a plots that visualizes the \eqn{-log_{10}(p)} QTL values (y-axis) 
  along the markers on the genome (x-axis), parameterized based on user-defined settings.
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
The \code{markers} should be ordered sequentially. 
The names of makers should be consistent over \code{markers} and \code{qtlProfiles}.
}
\seealso{ 
Use \code{\link{markers}} as example data set or use \code{loadData} to load your own data.
Use \code{\link{qtlMapTwoPart}} to calculate \code{qtlProfiles}.
Use \code{\link{qtlThreshold}} to estimate \code{qtlThres} QTL significance 
threshold .
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
                                             
##calculate the two part qtl
qtlProfiles <- qtlMapTwoPart(genotypes=genotypes, traits=traits[1:5,], spike=4)

##set the qtl threshold
qtlThres    <- 3.79

##OR: estimate the threshold yourself
#qtlThres   <- qtlThreshold(genotypes, traits, spike=4)
                                                                                                                       
##plot the qtlProfiles
qtlPlot(markers, qtlProfiles, qtlThres)
}