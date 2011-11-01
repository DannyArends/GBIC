\name{genotypes}
\docType{data}
\alias{genotypes}
\title{Example genotype data}
\description{
Example matrix of \code{genotypes} for each marker (rownames) and 162 Arabidopsis thaliana 
individuals (columnnames), as numeric values 1, 2 or NA when missing:
\preformatted{
> data(genotypes)
> genotypes[1:5,1:5]}
\tabular{lrrrrr}{
                      \tab X1   \tab X3   \tab X4     \tab X5     \tab X6   \cr
PVV4                    \tab 1    \tab 1      \tab 2      \tab 1    \tab 2    \cr
AXR-1                   \tab 1      \tab 1    \tab 2      \tab 1      \tab 2    \cr
HH.335C-Col             \tab 1      \tab 1    \tab 1      \tab 1      \tab 2    \cr
DF.162L/164C-Col        \tab 1      \tab 1    \tab 1      \tab 1      \tab 2    \cr
EC.480C               \tab 1    \tab 1    \tab 1    \tab 1    \tab 2    \cr
}
}
\usage{data(genotypes)}
\format{ 117 markers by 162 samples (RILS). }
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
Joost Keurentjes <Joost.Keurentjes@wur.nl>,
Morris Swertz <m.a.swertz@rug.nl>, 
Ritsert Jansen <r.c.jansen@rug.nl>
}
\seealso{ 
  Use \code{\link{MetaNetwork}} for automated analysis of this data set as part 
  a genetic analysis protocol on metabolites.
  Use \code{\link{qtlMapTwoPart}} for the calculation of \code{qtlProfiles} form \code{genotypes}.  
}
\examples{
  ##load the data
  data(genotypes)
  
  ##view part of the the data
  genotypes[1:5,1:5]
}
\keyword{datasets}