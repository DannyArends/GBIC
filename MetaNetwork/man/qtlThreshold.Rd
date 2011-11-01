\name{qtlThreshold}
\alias{qtlThreshold}
\title{ Estimate QTL significance threshold. }
\description{
  Simulation test to estimate empirical threshold for QTL significance. The trait 
  values are simulated using the median number of noise, median mean and 
  standard deviation from the observed trait data under the null hypothesis of no QTL. 
  At each simulation test, the highest absolute \eqn{-log_{10}(p)} value is recorded. 
  The threshold is set at desired alpha level (i.e. take the \eqn{-log_{10}(p)} value at the
  alpha proportion position of the permutation set).
}
\usage{
qtlThreshold(genotypes, traits, spike, 
             n.simulations = 1000, alpha = 0.05)
}
\arguments{
  \item{genotypes}{ matrix of genotypes for each marker (rownames) and individual
  (columnnames), as numeric values 1, 2 or NA when missing.
  See \code{\link{genotypes}} example data. }
  \item{traits}{ matrix of phenotypes for each trait (rownames) and individual 
  (columnnames), as numeric or NA when missing.
  See \code{\link{traits}} example data.}
  \item{spike}{ numeric cut-off value to separate absent (qualitative) from 
  available (quantitative) trait abundance. }
  \item{n.simulations}{ (optional) number of simulations. Default is 1000 times. }
  \item{alpha}{ (optional) numeric alpha level for the threshold. Default is 0.05. }
}
\value{
  Returns the \eqn{-log_{10}(p)} significance threshold value for QTLs.
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
The individual columns of \code{genotypes} and \code{traits} must have the same order.
The names of individuals should be consistent over \code{genotypes} and \code{traits}.
}
\seealso{
Use \code{\link{genotypes}} and \code{\link{traits}} as example data sets or use 
\code{\link{loadData}} to load your own data. 
Use \code{\link{qtlFDR}} to estimate QTL significance threshold 
to control false discovery rate FDR.
Use \code{\link{MetaNetwork}} for automated appliction of this function as part 
a genetic analysis protocol on metabolites.
}
\examples{
## load the example data provided with this package                          
data(genotypes)
data(traits)                                                       
                                             
##OR: load your own data                     
#genotypes  <- loadData("genotypes.csv")
#traits     <- loadData("traits.csv")  
                                             
##estimate qtl threshold for significance with low count of simulations (advised: 1000)
qtlThres    <- qtlThreshold(genotypes, traits, 4, n.simulations = 3)

##show the threshold
qtlThres
}