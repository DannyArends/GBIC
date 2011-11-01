\name{qtlFDR}
\alias{qtlFDR}
\title{ Estimate QTL thresholds for false discovery rate (FDR) }
\description{
  Calculate thresholds to control the false discovery rate in QTL analysis. 
  At given fdrThres(desired proportion of false positives), estimate the corresponding qtlThreshold (qtlThres). 
  And/or at given qtlThres (desired treshold for significance), estimate the corresponding proportion of false positives (fdrThres). 
}
\usage{
qtlFDR(qtlProfiles, fdrThres=0.05, qtlThres=NULL ) 
}
\arguments{
  \item{qtlProfiles}{ matrix of QTL mapping of traits (rownames) 
  to markers (columnnames), as \eqn{-log_{10}(p)} values.
  See \code{\link{qtlProfiles}} example data.}
  \item{fdrThres}{ the desired qvalue proportion of false positives incured (called the false discovery rate). Default is 0.05.
  See \code{\link{qvalue}} package.}
  \item{qtlThres}{ numeric \eqn{-log_{10}(p)} desired threshold value for significant QTLs.
  See \code{\link{qtlThreshold}} function.}
}
\value{
A matrix with three columns and two rows:
\item{c: qValue}{proportion of false positives (fdr).}
\item{c: pValue}{\eqn{10^{-1*qtlThres}} transformation of qtlThres/qtlProfiles to match \code{\link{qvalue}} package.}
\item{c: -log10P}{\eqn{-log_{10}(p)} estimation of \code{qtlThres}.}
\item{r: fdrThres}{a pValue and \eqn{-log_{10}(p)} qtlThres for the given \code{fdrThres}.}
\item{r: qtlThres}{a qValue (fdr) and pValue for the given \code{qtlThres}.}
Example:
\tabular{lrrr}{
     \tab      qValue \tab       pValue \tab -log10P  \cr
[1,] \tab 0.050000000 \tab 0.0809722708 \tab 1.091664 \cr
[2,] \tab 0.001344755 \tab 0.0006845554 \tab 3.164591 \cr
}
If \code{fdrThres} or \code{qtlThres} is NULL then the respective rows are ommitted.
}
\source{
  Storey, J. D. & Tibshirani, R. Statistical significance for genomewide 
  studies. Proc. Natl. Acd. Sci. USA 100, 9440-9445 (2003).
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
  Use \code{\link{qtlThreshold}} to estimate an QTL threshold based on simulation. 
  Use \code{\link{qtlMapTwoPart}} to calculate \code{\link{qtlProfiles}}.
  Use \code{\link{MetaNetwork}} for automated application of this function as part a genetic analysis protocol on metabolites.
}
\examples{
## load the example data provided with this package                         
data(markers)   
data(genotypes)
data(traits)                                                       
                                             
##OR: load your own data                     
#markers        <- loadData("markers.csv")
#genotypes      <- loadData("genotypes.csv")
#traits         <- loadData("traits.csv")  
                                             
##calculate the two part qtl
qtlProfiles     <- qtlMapTwoPart(genotypes=genotypes, traits=traits, spike=4)
  
##set the qtl threshold
qtlThres        <- 3.79

##OR: estimate the threshold yourself
#qtlThres       <- qtlThreshold(genotypes, traits, spike=4, n.simulations=3)

##estimate FDR
qtlFDR          <- qtlFDR(qtlProfiles, fdrThres=0.05, qtlThres=qtlThres)

##show FDR for both fdrThres of 0.05 and qtlThres of qtlThres
qtlFDR
}