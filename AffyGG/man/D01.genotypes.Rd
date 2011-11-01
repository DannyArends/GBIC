\name{genotypes}
\docType{data}
\alias{genotypes}
\title{Example genotype data}
\description{
\code{genotypes}: Example matrix of \code{genotypes} for each marker (rownames) and 30 mice  
individuals (columnnames), as numeric values 1, 2 or NA when missing, indicating inheritance
from the first or second individual of the (starting) parents:
\preformatted{
> data(genotypes)
> genotypes[1:5,1:5]}\cr
\tabular{lrrrrr}{
                    \tab BXD1 \tab BXD2 \tab BXD5   \tab BXD6   \tab BXD8   \cr
D01Msw003             \tab 1    \tab 1    \tab 2      \tab 2      \tab 2    \cr
D01Msw011             \tab 1    \tab 1    \tab 2      \tab 2      \tab 2    \cr
D1Mit294              \tab 1    \tab 1    \tab 2      \tab 2      \tab 2    \cr
D1Mit67               \tab 1    \tab 1    \tab 1      \tab 2      \tab 2    \cr
D1Mit430              \tab 1    \tab 1    \tab 1      \tab 2      \tab 2    \cr
}
}
\usage{
data(genotypes)
}
\format{ 
\code{genotypes}: 705 markers by 30 samples (RILs). 
}
\source{ 
  
  Alberts R, Terpstra P, Bystrykh LV, de Haan G, Jansen RC. A statistical
  multiprobe model for analyzing cis and trans genes in genetical genomics
  experiments with short-oligonucleotide arrays. Genetics. 2005 Nov;171(3):1437-9.
}

\author{ 
Rudi Alberts <r.alberts@rug.nl>,  
Ritsert Jansen <r.c.jansen@rug.nl>
}

\examples{
  ##load the data
  data(genotypes)
  
  ##view part of the the data
  genotypes[1:5,1:5]
}
\keyword{datasets}
