\name{markers}
\docType{data}
\alias{markers}
\title{Example marker data}
\description{
Example matrix of \code{markers} (rownames) for Arabidopsis thaliana and their chromosome numbers (column 1) 
and centi-Morgan (cM, column 2) positions, ordered by position:
\preformatted{
> data(markers)
> markers[1:5,]}
\tabular{lll}{  
                          \tab chr \tab cM   \cr
PVV4                \tab 1   \tab 0.0  \cr
AXR-1               \tab 1   \tab 6.4  \cr
HH.335C-Col         \tab 1   \tab 10.8 \cr
DF.162L/164C-Col    \tab 1   \tab 12.9 \cr
EC.480C             \tab 1   \tab 15.059 \cr
}
}
\usage{data(markers)}
\format{ 117 rows of markers and two columns for chromosome number and centi-Morgan position }
\source{ 
  Keurentjes JJB, FU J, de vos CHR, Lommen A, Hall RD, Bino RJ, van der Plas LHW, Jansen RC, Vreugdenhil D,
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
}
\examples{
  ##load the data
  data(markers)
  
  ##view part of the the data
  markers[1:5,]
}
\keyword{datasets}