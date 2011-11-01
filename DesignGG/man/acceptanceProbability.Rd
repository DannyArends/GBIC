\name{acceptanceProbability}
\alias{acceptanceProbability}
\title{ Compute the acceptance probability for each updated design }
\description{
 Compute the acceptance probability for each updated design. It depends on the
 current temperature value of simulated annealing process.
 This is a subfunction needed for \code{designGG}, but is not directly used.}
\usage{
    acceptanceProbability( designScore, newDesignScore, method, 
                           temperature )
}

\arguments{
  \item{designScore}{ 
      score of current design.
  }
  \item{newDesignScore}{ 
      score of updated design. 
  }
  \item{method}{ 
      either "SA" (simulated annealing) or "MH". 
  (Metropolis Hastings) 
  }
  \item{temperature}{ 
      current temperature in simulated annealing process. 
  }
}

\references{ 
  E. Wit and J. McClure. Statistics for Microarrays: Design, Analysis
      and Inference. (2004) Chichester: Wiley.  \cr
   Y. Li, R. Breitling and R.C. Jansen. Generalizing genetical 
     genomics: the added value from environmental perturbation, Trends Genet
     (2008) 24:518-524. \cr    
  Y. Li, M. Swertz, G. Vera, J. Fu, R. Breitling, and R.C. Jansen. designGG: 
  An R-package and Web tool for the optimal design of genetical genomics 
  experiments. (submitted) \cr
  http://gbic.biol.rug.nl/designGG  
}
\author{
  Yang Li <yang.li@rug.nl>, Gonzalo Vera <gonzalo.vera.rodriguez@gmail.com> \cr
  Rainer Breitling <r.breitling@rug.nl>, Ritsert Jansen <r.c.jansen@rug.nl> }
\seealso{ 
  \code{\link{designGG}}
}
\keyword{method}
