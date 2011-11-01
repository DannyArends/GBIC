\name{qtlProfiles}
\docType{data}
\alias{qtlProfiles}
\title{Example qtl profiles}
\description{
Example matrix of QTL mapping profiles of traits (rownames) 
  on to markers (columnnames), calculated from Arabidopsis thaliana \code{genotypes}
  and \code{traits} example datasets, 
as \eqn{-log_{10}(p)} values.
The +/- sign indicates the direction of the additive effect:                    
  values are positive if the QTL has higher metabolite abundance for individuals       
  carrying the genotype 2 than those carrying the genotype 1; values are negative 
  otherwise.
\preformatted{
>data(qtlProfiles)
>qtlProfiles[1:5,1:5]}
\tabular{lrrrrr}{
                       \tab  PVV4       \tab AXR-1        \tab HH.335C-Col \tab  DF.162L/164C-Col  \tab  EC.480C \cr
3-Hydroxypropyl        \tab  1.1003510  \tab 1.52622550   \tab 0.7823755   \tab  0.69234969        \tab  1.0440796 \cr
4-Hydroxybutyl         \tab  0.3007378  \tab 1.86596872   \tab 2.3859329   \tab  2.24693677        \tab  1.7367161 \cr
4-Methylsulfinylbutyl  \tab -0.1738533  \tab 0.08042058   \tab 0.1171418   \tab  0.02829792        \tab -0.1050181 \cr
3-Butenyl              \tab -0.1875411  \tab -1.16887454  \tab -0.8904340  \tab -0.80973239        \tab -0.9517352 \cr
3-Methylthiopropyl     \tab  0.5516505  \tab 0.33930888   \tab 0.3652804   \tab  0.36563859        \tab  0.2445772 \cr
}
This matrix can be produced by functions \code{\link{qtlMapTwoPart}} or \code{\link{MetaNetwork}}.
}
\usage{data(qtlProfiles)}
\format{ A matrix of QTL \eqn{-log_{10}(p)} values for 24 traits (rows) and 117 markers(columns). }
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
  Use \code{\link{qtlMapTwoPart}} to generate \code{qtlProfiles}.
  Use \code{\link{MetaNetwork}} for automated analysis of this data set as part a genetic analysis protocol on metabolites.
}
\keyword{datasets}