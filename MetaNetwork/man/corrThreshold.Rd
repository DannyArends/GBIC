\name{qtlCorrThreshold}
\alias{qtlCorrThreshold}
\title{ Estimate correlation threshold }
\description{
  Permutation test to estimate threshold for significant QTL correlation.
  In each permutation, the relation between genotype and trait are disturbed
  and then correlation is calculated. Thus, any detected genetic correlation is assumed to be false. 
  Compute the partial correlation coefficient in each permutation and record highest 
  value. Then after \code{n.permutations}, the threshold is set at desired \code{alpha} level, 
  Bonferroni corrected by the number of edges per trait (the number of traits-1).  
}
\usage{
qtlCorrThreshold(markers, genotypes, traits, spike, qtlThres, 
                 n.permutations = 10000, alpha = 0.05, method = "qtl")
}
\arguments{
  \item{markers}{ matrix of markers (rownames) and their chromosome numbers (column 1)
   and centi-Morgan positions (cM, column 2), ordered by position.
   See \code{\link{markers}} example data. }
  \item{genotypes}{ matrix of genotypes for each marker (rownames) and individual 
  (columnnames), as numeric values 1, 2 or NA when missing.
  See \code{\link{genotypes}} example data.}
  \item{traits}{ matrix of phenotypes for each trait (rownames) and individual 
  (columnnames), as numeric or NA when missing.
  See \code{\link{traits}} example data.}
  \item{spike}{ numeric cut-off value to separate absent (qualitative) from 
  available (quantitative) trait abundance. }
  \item{qtlThres}{ numeric \eqn{-log_{10}(p)} threshold value for significant QTLs. }
  \item{n.permutations}{  (optional) number of iterations of permutation. 
  Default is 10000 times. }
  \item{alpha}{ (optional) Bonferonni corrected significant level. Default is 0.05. }
  \item{method}{ (optional) correlation method string, "qtl" or "abundance"  
  for correlation on QTL profiles or 
  metabolite abundance level, respectively. 
  If \code{method} is "qtl" (default), 
  qtlCorrThreshold will call \code{\link{qtlCorrZeroOrder}} to calculate the correlation between 
  QTL profiles. Otherwise, when \code{method} is "abundance", qtlCorrThreshold will 
  use Spearman correlation with the \code{\link{cor}} function to calculate the 
  correlation between metabolite 
  abundance profiles.}
}
\value{
 Returns: threshold value for correlation coefficient.
 \item{corrPermutations}{a vector of the permutations of maximum, absolute correlation values.}
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
The \code{markers} should be ordered sequentially.
The names of markers, traits and individuals should be consistent over \code{markers}, 
\code{genotypes}, and \code{traits}.
}
\seealso{ 
Use \code{\link{markers}}, \code{\link{genotypes}} and \code{\link{traits}} as example data sets or use \code{\link{loadData}} to load your own data. 
Use \code{\link{qtlThreshold}} to estimate \code{qtlThres} QTL significance threshold.
Use \code{\link{qtlCorrZeroOrder}} and \code{\link{qtlCorrSecondOrder}} for QTL correlations that can be test against the resulting correlation threshold.
Use \code{\link{MetaNetwork}} for automated appliction of this function as part a genetic analysis protocol on metabolites.
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
                                             
##set the qtl threshold
qtlThres    <- 3.79

##OR: estimate the threshold yourself
#qtlThres   <- qtlThreshold(genotypes, traits, spike=4)

##estimate qtlCorrThreshold with low number of permutations (advised: 10000)
cat("warning: simulating correlations takes a few minutes")
flush.console()
\dontrun{
corrThres   <- qtlCorrThreshold(markers, genotypes, traits, spike=4, 
                                qtlThres=qtlThres, n.permutations=3) 

##show threshold
corrThres
}
}