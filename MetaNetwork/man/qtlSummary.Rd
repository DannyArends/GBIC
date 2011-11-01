\name{qtlSummary}
\alias{qtlSummary}
\title{ Summarize QTL effects }
\description{
  Summarize QTL profiles including peak position (chromosome, centi-Morgan), 
  support interval (calculated via the function \code{qtlSupportInterval}), 
  proportion of QTL variation explained by part 1 and 2 of the \code{\link{qtlMapTwoPart}}
  model, and allele substitution effect.
}
\usage{
  qtlSummary(markers, genotypes, traits, qtlProfiles, spike, qtlThres, 
             interval.dropoff = 1.5, filename = NULL)
}
\arguments{
  \item{markers}{ matrix of markers (rownames) and their chromosome numbers 
  (column 1) and centi-Morgan positions (cM, column 2), ordered by position.
  See \code{\link{markers}} example data. }
  \item{genotypes}{ matrix of genotypes for each marker (rownames) and 
  individual (columnnames), as numeric values 1, 2 or NA when missing.
  See \code{\link{genotypes}} example data.}
  \item{traits}{ matrix of phenotypes for each trait (rownames) and individual 
  (columnnames), as numeric or NA when missing.
  See \code{\link{traits}} example data. }
  \item{qtlProfiles}{ matrix of QTL mapping of traits (rownames) 
  to markers (columnnames), as \eqn{-log_{10}(p)} values.
  See \code{\link{qtlProfiles}} example data.}
  \item{spike}{ numeric cut-off value to separate absent (qualitative) from 
  available (quantitative) trait abundance.  }
  \item{qtlThres}{ numeric \eqn{-log_{10}(p)} threshold value for significant QTLs. }
  \item{interval.dropoff}{ (optional) drop-off value for QTL support 
  intervals. Default is 1.5. }
  \item{filename}{ (optional) path of the file where the qtlSummary is to be 
  stored. Default is NULL.}
}
\value{
  Returns a data frame with a QTL summary which contains the following headers:
  \item{traitName}{ name of trait.}   
  \item{QTLchr}{ the chromosome number where a QTL locates.}    
  \item{QTLmk}{ the marker where the trait maps to.}
  \item{QTLleftcm}{ the cM position of left border of the QTL support interval.}
  \item{QTLpeakcm}{ the cM position of the QTL peak.}    
  \item{QTLrightcm}{ the cM of right border of the QTL support interval.}
  \item{logp}{ the \eqn{-log_{10}(p)} value of a QTL.}
  \item{VarP1}{ the percentage of qualitative variance explained by a QTL.}
  \item{VarP2}{ the percentage of quantitative variance explained by a QTL.}
  \item{additive}{ the allele substitution effect (=half the difference of 
  metabolite abundance between genotype 1 and 2).} 
  See \code{\link{qtlSumm}} example data.  
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
The individual columns of genotypes and traits must have the same order. 
The markers should be ordered sequentially.
The names of markers, traits and individuals should be consistent over \code{markers}, 
\code{genotypes}, \code{traits} and \code{qtlProfiles}.
}
\seealso{ 
Use \code{\link{markers}}, \code{\link{genotypes}} and \code{\link{traits}} 
as example data sets or use \code{\link{loadData}} to load your own data. 
Use \code{\link{qtlMapTwoPart}} to calculate \code{qtlProfiles}.
Use \code{\link{qtlThreshold}} to estimate \code{qtlThres} QTL threshold for 
significance.
Use \code{\link{MetaNetwork}} for automated application of this function as 
part of a genetic analysis protocol on metabolites.
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
qtlProfiles <- qtlMapTwoPart(genotypes=genotypes, traits=traits, spike=4)

##set the qtl threshold
qtlThres    <- 3.79

##OR: estimate the threshold yourself
#qtlThres   <- qtlThreshold(genotypes, traits, spike=4)
                                                                                                                       
##summarize the qtlProfiles
qtlSumm     <- qtlSummary(markers, genotypes, traits, qtlProfiles, spike=4, 
                          qtlThres=qtlThres)

##show the summary
qtlSumm[1:5,]
}