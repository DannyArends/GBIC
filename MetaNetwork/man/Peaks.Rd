\name{peaks}
\docType{data}
\alias{peaks}
\alias{peaks2}
\title{Example peak data}
\description{Example matrix of mass/charge peaks (column1) for each trait (rownames)
that accompanies the untargeted \code{traits2} metabolite
abundance dataset from Arabidopsis thaliana:
\preformatted{
> data(peaks2)
> peaks2[1:5,]}
\tabular{lr}{
            \tab MassOverCharge \cr
LCavg.1539  \tab 378            \cr       
LCavg.1549  \tab 379            \cr      
LCavg.1555  \tab 753            \cr 
LCavg.1570  \tab 424            \cr
LCavg.1610  \tab 436            \cr      
}
}
\usage{data(peaks2)}
\format{ 20 rows of peaks with mass over charge values (column 1) for example data set \code{\link{traits2}}.
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
  Use \code{\link{MetaNetwork}} for automated analysis of this data set as part a genetic analysis protocol on metabolites.
  Use \code{\link{findPeakMultiplicity}} to relate multiple peaks per metabolite.
}
\examples{
## load the example data provided with this package
data(peaks2)

##view part of the the data
peaks2[1:5,]        
}
\keyword{datasets}