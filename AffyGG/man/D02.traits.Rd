\name{traits}
\docType{data}
\alias{traits}
\alias{traits2}
\title{Example traits data}
\description{
Example matrix of phenotypes for each trait (rownames) and sample/individual 
(columnnames), as numeric or NA when missing.

The \code{traits} dataset contains abundance of 16 probe
(gene expression) for 30 mice individuals. Probe set named  
\eqn{96243\_f\_at}{96243_f_at}.

\preformatted{
> data(traits)
> traits[1:5,1:5]}\cr
\tabular{lrrrrr}{

              \tab BXD01  \tab BXD02  \tab BXD05  \tab BXD06   \tab BXD08 \cr
96243\_f\_at1  \tab 8.8896 \tab 8.3223 \tab 9.7237 \tab 8.9042  \tab 9.4805 \cr
96243\_f\_at2  \tab 8.4075 \tab 7.7312 \tab 9.1945 \tab 8.0864  \tab 8.4424 \cr
96243\_f\_at3  \tab 7.9805 \tab 7.5844 \tab 8.6489 \tab 7.8963  \tab 8.1320 \cr
96243\_f\_at4  \tab 9.2778 \tab 8.6701 \tab 6.1207 \tab 8.6260  \tab 8.9784 \cr
96243\_f\_at5  \tab 10.3764 \tab 9.7419 \tab 7.8298 \tab 9.9070 \tab 10.4877 \cr
}

The \code{traits2} dataset contains abundance of 16 probe
(gene expression) for 30 mice individuals. Probe set named  
\eqn{160371\_at}{160371_at}. \cr

\preformatted{
> data(traits2)
> traits2[1:5,1:5]}\cr
\tabular{lrrrrr}{

            \tab   BXD01 \tab   BXD02 \tab   BXD05 \tab   BXD06 \tab   BXD08 \cr
160371\_at1 \tab 12.2860 \tab 12.2125 \tab 12.3105 \tab 11.9035 \tab 12.4335 \cr
160371\_at2 \tab  4.6550 \tab  4.7245 \tab  4.6900 \tab  4.7735 \tab  4.7130 \cr
160371\_at3 \tab  6.5155 \tab  6.4180 \tab  6.4195 \tab  7.1330 \tab  7.0515 \cr
160371\_at4 \tab 12.3980 \tab 12.4715 \tab 12.4520 \tab 12.0920 \tab 12.1355 \cr
160371\_at5 \tab  6.3675 \tab  6.6170 \tab  5.9580 \tab  5.5850 \tab  5.8235 \cr
}

}
\usage{
data(traits)
data(traits2)
}
\format{ The number of traits (rows) by the number of samples (columns, RILS). }

\examples{
  ##load the identified data set
  data(traits)        
  ##view part of the the data
  traits[1:5,1:5]

  ##load the identified data set
  data(traits2)        
  ##view part of the the data
  traits2[1:5,1:5]
}
\source{   
  Alberts R, Terpstra P, Bystrykh LV, de Haan G, Jansen RC. A statistical
  multiprobe model for analyzing cis and trans genes in genetical genomics
  experiments with short-oligonucleotide arrays. Genetics. 2005 Nov;171(3):1437-9. Epub 2005 Aug 3. 
}

\author{
Rudi Alberts <r.alberts@rug.nl>,
Ritsert Jansen <r.c.jansen@rug.nl>
}
\keyword{
datasets
}
