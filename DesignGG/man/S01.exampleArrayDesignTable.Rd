\name{exampleArrayDesignTable}
\docType{data}
\alias{exampleArrayDesignTable}
\title{Example output of ArrayDesignTable data}
\description{
\code{exampleArrayDesignTable}: Example data of \code{exampleArrayDesignTable} for a 
hypothetical dual-channel microarray experiment in which there are 100 recombinant
inbred lines (RIL) and 27 arrays available. Two environmental factors 
(temperature and cell type) are considered in this experiment. There are three 
levels for temperature (15, 24 and 29) and four levels for cell types (A,B,C,D).
This table tells how to pair samples into arrays.  

\preformatted{
> data(exampleArrayDesignTable)
> exampleArrayDesignTable[1:5,]}\cr
\tabular{lrr}{

            \tab Channel1 \tab Channel2 \cr 
array1      \tab RIL28    \tab RIL92 \cr    
array2      \tab RIL70    \tab RIL47 \cr    
array3      \tab RIL22    \tab RIL89 \cr    
array4      \tab RIL45    \tab RIL15 \cr    
array5      \tab RIL52    \tab RIL41 \cr
}
}
\usage{
data(exampleArrayDesignTable)
}
\format{ 
\code{exampleArrayDesignTable}: 27 arrays by two channels. 
}

\author{ 
 Yang Li <yang.li@rug.nl> \cr
 Rainer Breitling <r.breitling@rug.nl> \cr
 Ritsert Jansen <r.c.jansen@rug.nl>
}

\examples{
  ##load the data
  data(exampleArrayDesignTable)
  
  ##view part of the the data
  exampleArrayDesignTable[1:5,]
}
\keyword{datasets}
