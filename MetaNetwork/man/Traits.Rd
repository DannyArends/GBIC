\name{traits}
\docType{data}
\alias{traits}
\alias{traits2}
\title{Example traits data}
\description{
The \code{traits} dataset contains measured abundance of 18 glucosinates and 6 flavonols
for 162 Arabidopsis thaliana individuals. The \code{traits2} dataset contains abundance of 20 unidentified 
mass peaks (mass peak numbers) for 162 Arabidopsis thaliana individuals.
Example matrice of phenotypes for each trait (rownames) and sample/individual 
(columnnames), as numeric or NA when missing:
\preformatted{
> data(traits)
> traits[1:5,1:5]}
\tabular{lrrrrr}{
                       \tab X1  \tab X3  \tab X4   \tab X5    \tab X6    \cr
3-Hydroxypropyl        \tab NA  \tab 942 \tab 2402 \tab 602   \tab 213   \cr
4-Hydroxybutyl         \tab NA  \tab 4   \tab 10   \tab 183   \tab 198   \cr
4-Methylsulfinylbutyl  \tab NA  \tab 55  \tab 62   \tab 13386 \tab 1671  \cr
3-Butenyl              \tab NA  \tab 84  \tab 32   \tab 18    \tab 4339  \cr
3-Methylthiopropyl     \tab NA  \tab 3108\tab 569  \tab 4     \tab 7     \cr
}
}
\usage{data(traits)
data(traits2)}
\format{ the number of traits (rows) by the number of samples (columns, RILS). }
\seealso{ 
  Use \code{\link{MetaNetwork}} for automated analysis of this data set as part 
  a genetic analysis protocol on metabolites.
  Use \code{\link{qtlMapTwoPart}} for the calculation of \code{qtlProfiles} form \code{traits}.  
}
\examples{
  ##load the identified data set
  data(traits)
  
  ##view part of the the data
  traits[1:5,1:5]
  
  ##load the unidentified data set
  data(traits2)
  
  ##view part of the the data
  traits2[1:5,1:5]  
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
Joost Keurentjes <Joost.Keurentjes@wur.nl>, 
Morris Swertz <m.a.swertz@rug.nl>, 
Ritsert Jansen <r.c.jansen@rug.nl>
}
\keyword{datasets}