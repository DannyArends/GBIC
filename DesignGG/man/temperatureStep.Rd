\name{temperatureStep}
\alias{temperatureStep}
\title{Calculate the temperature decreasing step for simulated annealing process}

\description{   
    Calculate the temperature decreasing step for simulated annealing process.
    This is a subfunction needed for \code{designGG}, but is not directly used.
}
\usage{
temperatureStep(startTemp, maxTempStep, endTemp, nIterations)
}
\arguments{
  \item{startTemp}{
      starting temperature of simulated annealing process.
  }
  \item{maxTempStep}{ 
      maximum  temperature decreasing step for simulated 
      annealing process. The parameter ensures that the multiplicative cooling 
      factor is not smaller than this value. If nIterations is too small, the preferred 
      final temperature (endTemp) may not be reached. See Wit and McClure (2004) 
      for details. 
  }
  \item{endTemp}{ 
      ending temperature of simulated annealing process.
      An important optimization parameter. Setting this parameter closer to zero.
      See Wit and McClure (2004) for details 
  }
  \item{nIterations}{ 
      number of iterations in the simulated annealing method.
  }
}
\value{
   A temperature decreasing step in the simulated annealing process.
 
} 
\references{
  Y. Li, M. Swertz, G. Vera, J. Fu, R. Breitling, and R.C. Jansen. designGG: 
  An R-package and Web tool for the optimal design of genetical genomics 
  experiments. (submitted) \cr
  http://gbic.biol.rug.nl/designGG \cr
  Y. Li, R. Breitling and R.C. Jansen. Generalizing genetical 
     genomics: the added value from environmental perturbation, Trends Genet
     (2008) 24:518-524. \cr
  E. Wit and J. McClure. Statistics for Microarrays: Design, Analysis
      and Inference. (2004) Chichester: Wiley.    
}
\author{  
  Yang Li <yang.li@rug.nl>, Gonzalo Vera <gonzalo.vera.rodriguez@gmail.com> \cr
  Rainer Breitling <r.breitling@rug.nl>, Ritsert Jansen <r.c.jansen@rug.nl> 
}
\seealso{ 
  \code{\link{designGG}}
}
\keyword{method}
