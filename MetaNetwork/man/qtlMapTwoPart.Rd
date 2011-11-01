\name{qtlMapTwoPart}
\alias{qtlMapTwoPart}
\title{ QTL test statistic using two-part model }
\description{
  QTL model for a mixed phenotype with qualitative and quantitative variance.
  
  Part one tests whether the presence/absence of metabolites has a genetic 
  basis, i.e. whether different genotype classes at a given marker differ in 
  their number of non-zero observations. 
  
  Part two tests whether 
  quantitative variation in non-zero abundance has a genetic basis: whether 
  the non-zero observations for each of these genotype classes at a given 
  marker differ in mean abundance. 
  
  The 'p-value' of the QTL is computed 
  as the product of the two 'p-values' in the two parts. 
}
\usage{
qtlMapTwoPart(genotypes, traits, spike, filename=NULL) 
}
\arguments{
  \item{genotypes}{ matrix of genotypes for each marker (rownames) and individual 
  (columnnames), as numeric values 1, 2 or NA when missing.
  See \code{\link{genotypes}} example data.}
  \item{traits}{ matrix of phenotypes for each trait (rownames) and individual 
  (columnnames), as numeric or NA when missing.
  See \code{\link{traits}} or \code{\link{traits2}} example data. }
  \item{spike}{ numeric cut-off value to separate absent (qualitative) from available (quantitative) trait abundance. }
  \item{filename}{ (optional) path of the file where the QTLs are to be stored. Default is NULL. }
}
\details{
  The analysis model for \code{traits} data with both qualitative and quantitative variance. 
  If the distribution of \code{traits} is mixed, a two part model is used. 
  
  The first part describes a binominal model which tests 
  for association of markers with observed noise value. For each trait: 
  Let \eqn{y_i} denote the mass intensity for \eqn{i^{th}} RIL and let \eqn{z_i=0} if \eqn{y_i \leq} \code{spike} and \eqn{z_i=1} if \eqn{y_i >}\code{spike}.
  We then test each marker for significant differences between the two genotypes 
  for the probability of presence of the trait (mass) peak: 
  \eqn{ H_0: p \lbrace z=1\arrowvert g=1 \rbrace = p \lbrace z=1 \arrowvert g=2 \rbrace}
  versus the alternative hypothesis
  \eqn{ H_1: p \lbrace z=1 \arrowvert g=1 \rbrace \neq p \lbrace z=1 \arrowvert g=2 \rbrace},
  where \eqn{g} is the genotype (1 or 2) of a marker under analysis. 
  
  The second part describes a parametric model which tests for association of 
  markers with quantitative variances for those individuals where \eqn{y_i >}\code{spike}.
  Under the assumption of normal distribution, we test each marker for 
  significant differences in the mean values between two genotypes: 
  \eqn{ H_0: u \lbrace g=1 \rbrace = u \lbrace g=2 \rbrace }
  versus the alternative hypothesis
  \eqn{ H_1: u \lbrace g=1 \neq u \lbrace g=2 \rbrace }. 
  
  When there is binary data only (no quantitative data) or quantitative data only (no spike), the 'p-value' of the missing part is set to one. 
  The p value of the two-part model is then determined by the multiple of the p 
  values from the two parts.
}
\value{
  Returns a matrix of QTL mapping of traits (rownames) 
  to markers (columnnames), as \eqn{-log_{10}(p)} values. 
  The \eqn{-log_{10}(p)} values are assigned negative if traits for individuals carrying genotype "1" 
  are larger on average than traits for individuals carrying genotype "2". 
  Otherwise, the \eqn{-log_{10}(p)} value is positive.
}
\source{
  Broman, K.W. Mapping quantitative trait loci in the case of a spike in the phenotype 
  distribution. Genetics 163, 1169-1175 (2003). 
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
The individual columns of \code{genotypes} and \code{traits} must have the same order.
The names of individuals must be consistent over \code{genotypes} and \code{traits}.
}
\seealso{ 
  Use \code{\link{genotypes}} and \code{\link{traits}} as example data sets or 
  use \code{\link{loadData}} to load your own data.
  Use \code{\link{MetaNetwork}} for automated application of this function as 
  part a genetic analysis protocol on metabolites.
}
\examples{
## load the example data provided with this package
data(genotypes)
data(traits)

##OR: load your own data
#genotypes  <- loadData("genotypes.csv")
#traits     <- loadData("traits.csv")

##calculate the two part qtl using a cutoff spike of 4
qtlProfiles <- qtlMapTwoPart(genotypes=genotypes, traits=traits[1:5,], spike=4) 

#show part of the results
qtlProfiles[1:5,1:5]
}