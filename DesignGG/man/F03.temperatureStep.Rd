\name{temperatureStep}
\alias{temperatureStep}
\title{Calculate the temperature decreasing step for simulated annealing process}

\description{   
    Calculate the temperature decreasing step for simulated annealing process.
    This is a subfunctions needed for g3d, but is not directly used.
}
\usage{
temperatureStep(startTemp, maxTempStep, endTemp, nIterations)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{startTemp}{
      starting temperature of simulated annealing process.
  }
  \item{maxTempStep}{ 
      maximum  temperature decreasing step for simulated 
      annealing process. The parameter ensures that the multiplicative cooling 
      factor is not smaller than that. If nIterations is too small, the preferred 
      final temperature (endTemp) may not be reached.See Wit and McClure (2004) 
      for details 
  }
  \item{endTemp}{ 
      ending temperature of simulated annealing process.
      An important optimization parameter. Setting this parameter closer to zero
      See Wit and McClure (2004) for details 
  }
  \item{nIterations}{ 
      Number of iterations of the simulated annealing method
  }
}
\references{
    Y. Li, R. Breitling and R.C. Jansen "Generalizing genetical 
     genomics: the added value from environmental perturbation", Trends Genet
     (2008) 24:518-524.
    E. Wit and J. McClure (2004) "Statistics for Microarrays: Design, Analysis
      and Inference" Chichester: Wiley. 
}
\author{  
 Yang Li <yang.li@rug.nl>, Gonzalo Vera <g.vera.rodriguez@rug.nl> \cr
 Rainer Breitling <r.breitling@rug.nl>, Ritsert Jansen <r.c.jansen@rug.nl> 
}
\keyword{method}
