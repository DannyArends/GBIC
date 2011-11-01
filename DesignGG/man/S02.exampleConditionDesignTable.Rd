                                                                                                            \name{exampleArrayDesignTable}
\docType{data}
\alias{exampleConditionDesignTable}
\title{Example ConditionDesignTable data}
\description{
\code{exampleConditionDesignTable}: 
Example data of \code{exampleConditionDesignTable} for a hypothetical 
dual-channel microarray experiment in which there are 100 recombinant 
inbred lines (RIL) and 27 arrays available. Two environmental factors 
(temperature and cell type) are considered in this experiment. There are three 
levels for temperature (15,24 and 29) and four levels for cell types (A,B,C,D).
This table tells how to allocate samples into 12 (=3*4) different conditions. 
On average there are 4.5 (27*2/12) samples per condition. 

\preformatted{
> data(exampleConditionDesignTable)
> exampleConditionDesignTable[1:5,]}\cr
\tabular{lrrrrrrr}{
            \tab Temperature  \tab Cell Type   \tab \tab \tab \tab \tab   \cr 
condition1  \tab 15 \tab A \tab RIL28 \tab RIL81 \tab RIL18 \tab RIL61 \tab \cr  
condition2  \tab 24 \tab A \tab RIL72 \tab RIL40 \tab RIL83 \tab RIL44 \tab RIL10 \cr 
condition3  \tab 29 \tab A \tab RIL22 \tab RIL89 \tab RIL3 \tab RIL30 \tab RIL58 \cr 
condition4  \tab 15 \tab B \tab RIL70 \tab RIL47 \tab RIL4 \tab RIL59 \tab \cr  
condition5  \tab 24 \tab B \tab RIL93 \tab RIL97 \tab RIL49 \tab RIL14 \tab \cr
}
}
\usage{
data(exampleConditionDesignTable)
}
\format{ 
\code{exampleConditionDesignTable}: 12 combination of conditions from three 
temepratures and four cell types.  
}

\author{ 
 Yang Li <yang.li@rug.nl> \cr
 Rainer Breitling <r.breitling@rug.nl> \cr
 Ritsert Jansen <r.c.jansen@rug.nl>
}

\examples{
  ##load the data
  data(exampleConditionDesignTable)
  
  ##view part of the the data
  exampleConditionDesignTable[1:5,]
}
\keyword{datasets}
