\name{corrZeroOrder and corrSecondOrder}
\docType{data}
\alias{corrZeroOrder}
\alias{corrSecondOrder}
\title{Example correlation data}
\description{
Example of a QTL correlation matrix of \code{\link{traits}} that can be 
produced using functions \code{\link{qtlCorrZeroOrder}} or \code{\link{qtlCorrSecondOrder}}.
\preformatted{
> data(corrZeroOrder)
> corrZeroOrder[1:3,1:3]}
\tabular{lrrr}{
                       \tab 3-Hydroxypropyl \tab 4-Hydroxybutyl \tab 4-Methylsulfinylbutyl \cr
3-Hydroxypropyl        \tab       1.0000000 \tab      0.8205883 \tab             0.8471439 \cr
4-Hydroxybutyl         \tab       0.8205883 \tab      1.0000000 \tab             0.9182831 \cr
4-Methylsulfinylbutyl  \tab       0.8471439 \tab      0.9182831 \tab             1.0000000 \cr    
}
}
\usage{data(corrZeroOrder)
data(corrSecondOrder)}
\format{ 24x24 matrix of correlation coefficients between \code{\link{traits}}.
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
  Use \code{\link{qtlCorrZeroOrder}},  \code{\link{qtlCorrSecondOrder}} or \code{\link{cor}} to reproduce this data set.
  Use \code{\link{MetaNetwork}} to produce this dataset as part a genetic analysis protocol on metabolites.
}
\examples{
## load the example data provided with this package
data(corrZeroOrder)
data(corrSecondOrder)
 
##OR: generate your own correlation data using MetaNetwork  
#genotypes  <- loadData("genotypes.csv") 
#traits     <- loadData("traits.csv") 
#markers    <- loadData("markers.csv")  
#
#MetaNetwork (markers=markers, genotypes=genotypes, traits=traits, spike=4) 
             
##view the correlation data
corrZeroOrder
corrSecondOrder              
}
\keyword{datasets}