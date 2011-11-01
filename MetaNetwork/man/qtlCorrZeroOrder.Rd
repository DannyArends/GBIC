\name{qtlCorrZeroOrder}
\alias{qtlCorrZeroOrder}
\title{Calculate the zero-order correlation}
\description{Calculate the zero-order correlation on QTL profiles.}
\usage{qtlCorrZeroOrder(markers, qtlProfiles, qtlThres, filename=NULL)}
\arguments{
  \item{markers}{ matrix of markers (rownames) and their chromosome numbers (column 1) 
  and centi-Morgan positions (cM, column 2), ordered by position.
  See \code{\link{markers}} example data. }
  \item{qtlProfiles}{ matrix of QTL mapping of traits (rownames) 
  to markers (columnnames), 
  as \eqn{-log_{10}(p)} values.
  See \code{\link{qtlProfiles}} example data.}
  \item{qtlThres}{ numeric \eqn{-log_{10}(p)} threshold value for significant QTLs. }
  \item{filename}{ (optional) path of the file where the correlations are to be stored. Default NULL. }
}
\details{
  QTL support intervals are determined (via \code{qtlSupportInterval} with \code{interval.dropoff = 1.5}) 
  and the \eqn{-log_{10}(p)} values outside of the 
  borders of these intervals are set to zero. Pairwise correlation coefficients between
  any two traits are then calculated as  
  
  \deqn{r_{xy} = \frac{2\displaystyle\sum_{i=1}^n x_i*y_i}{\displaystyle\sum_{i=1}^n x_i^2+\displaystyle\sum_{i=1}^n y_i^2} }
  
  where \eqn{r_{xy}} is the correlation coefficient between qtlProfiles \eqn{x} 
  and \eqn{y} and \eqn{i (i=1...n)} is a marker. \eqn{x_i} and \eqn{y_i} 
  represent \eqn{-log_{10}(p)} QTL profile values for marker \eqn{i}.
}
\value{
  Returns a matrix of correlation coefficients.
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
\note{
The markers should be ordered sequentially. 
The names of markers and traits should be consistent over \code{qtlProfiles} and \code{markers}.
}
\seealso{ 
Use \code{\link{markers}} as example data set or use \code{\link{loadData}} to load your own data.
Use \code{\link{qtlSupportInterval}} to calculate support intervals.
Use \code{\link{qtlMapTwoPart}} to calculate \code{qtlProfiles}.
Use \code{\link{qtlThreshold}} to estimate \code{qtlThres} QTL significance threshold .
Use \code{\link{MetaNetwork}} for automated application of this function as part a genetic analysis protocol on metabolites.
}
\examples{
## load the example data provided with this package                         
data(markers)   
data(genotypes)
data(traits)                                                       
                                             
##OR: load your own data                     
#markers     <- loadData("markers.csv")
#genotypes   <- loadData("genotypes.csv")
#traits      <- loadData("traits.csv")  
                                             
##calculate the two part qtl
qtlProfiles  <- qtlMapTwoPart(genotypes=genotypes, traits=traits[1:5,], spike=4)
  
##set the qtl threshold
qtlThres     <- 3.79

##OR: estimate the threshold yourself
#qtlThres    <- qtlThreshold(genotypes, traits, spike=4)
  
##calculate zero order correlation
qtlZeroOrder <- qtlCorrZeroOrder(markers, qtlProfiles, qtlThres)

##show the correlations
qtlZeroOrder[1:5,1:5]
}