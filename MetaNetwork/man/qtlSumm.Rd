\name{qtlSumm}
\docType{data}
\alias{qtlSumm}
\title{Example QTL summary report}
\description{
Example of the QTL summary for each trait as data frame that can be produced 
by functions \code{\link{qtlSummary}} or \code{\link{MetaNetwork}}.
\preformatted{
> data(qtlSumm)
> qtlSumm[1:5,]}
\tabular{llrrrrrr}{
   \tab traitName             \tab QTLchr \tab QTLmk   \tab QTLleftcm \tab QTLpeakcm \tab QTLrightcm \tab logP\cr
1  \tab 3-Hydroxypropyl       \tab 4      \tab GA1     \tab 3.6640  \tab 9.027   \tab 11.1280  \tab 11.1 \cr
2  \tab 3-Hydroxypropyl       \tab 5      \tab GH.117C \tab 32.4675 \tab 35.356  \tab 37.6390  \tab 13.5 \cr
3  \tab 4-Hydroxybutyl        \tab 3      \tab DF.77C  \tab 0.0000 \tab 0.000   \tab 9.4365  \tab  3.9 \cr
4  \tab 4-Hydroxybutyl        \tab 4      \tab GA1     \tab 3.6640  \tab 9.027   \tab 16.5810  \tab 4.5 \cr
5  \tab 4-Methylsulfinylbutyl \tab 5      \tab GH.117C \tab 32.356  \tab 35.356  \tab 37/6390  \tab  16.3 
}
\tabular{lrrr}{
(cont'd)   \tab VarP1 \tab VarP2 \tab additive \cr
(1 cont'd) \tab 0.0 \tab  26.0 \tab  -2784.0 \cr
(2 cont'd) \tab   0.0 \tab  31.2 \tab  -3088.8 \cr
(3 cont'd) \tab   10.4 \tab   1.3 \tab     -16.9 \cr
(4 cont'd) \tab 5.5 \tab  6.6 \tab     28.0 \cr
(5 cont'd) \tab   12.0 \tab  34.2 \tab  61.3
}

\code{qtlSumm} has the following headers:
\itemize{
    \item{\code{traitName}}{ name of trait.}   
    \item{\code{QTLchr}}{ the chromosome number where a QTL locates.}    
    \item{\code{QTLmk}}{ the marker where the trait maps to.}
    \item{\code{QTLleftcm}}{ the cM position of left border of the QTL support interval.}
    \item{\code{QTLpeak}}{ the cM position of the QTL peak.}    
    \item{\code{QTLrightcm}}{ the cM of right border of the QTL support interval.}
    \item{\code{logp}}{ the \eqn{-log_{10}(p)} value of a QTL.}
    \item{\code{VarP1}}{ the percentage of qualitative variance explained by a QTL.}
    \item{\code{VarP2}}{ the percentage of quantitative variance explained by a QTL.}
    \item{\code{additive}}{ the allele substitution effect (=half the difference of metabolite abundance between genotype 1 and 2).}}
}

\usage{data(qtlSumm)}
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
  Use \code{\link{qtlSummary}} to generate \code{qtlSumm}.
  Use \code{\link{MetaNetwork}} to produce this dataset as part a genetic 
  analysis protocol on metabolites.
}
\examples{
## load the example data provided with this package
data(qtlSumm)

##OR: generate your own qtlSummary using MetaNetwork
#genotypes  <- loadData("genotypes.csv") 
#traits     <- loadData("traits.csv") 
#markers    <- loadData("markers.csv")  
#
#MetaNetwork (markers=markers, genotypes=genotypes, traits=traits[1:2,], spike=4) 
             
##view the summary
qtlSumm              
}
\keyword{datasets}