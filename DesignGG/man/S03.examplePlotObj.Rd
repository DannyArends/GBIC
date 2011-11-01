\name{exampleArrayDesignTable}
\docType{data}
\alias{examplePlotObj}
\title{Example PlotObj data}
\description{
\code{examplePlotObj}: 
Example data of \code{examplePlotObj} for plot all scores and cooling at each 
iteration during simulated annealing process.

\preformatted{
> data(examplePlotObj)
> exampleConditionDesignTable
> plotAllScores(examplePlotObj)}\cr
}
\usage{
data(examplePlotObj)
}
\format{ 
\code{examplePlotObj}: a list which contains the following elements:
(1)scores (2)cooling (3) startTemp  (4)temprature 
(5)temperature.step (6)nIterations (7)optimality.  
}

\author{ 
 Yang Li <yang.li@rug.nl> \cr
 Rainer Breitling <r.breitling@rug.nl> \cr
 Ritsert Jansen <r.c.jansen@rug.nl>
}

\examples{
  ##load the data
  data(examplePlotObj)
  
  ##make plot based on the data
  plotAllScores(examplePlotObj)
  
}
\keyword{datasets}
