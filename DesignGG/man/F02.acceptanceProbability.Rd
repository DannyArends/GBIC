\name{acceptanceProbability}
\alias{acceptanceProbability}
\title{ Compute the acceptance probability for each updated design }
\description{
 Compute the acceptance probability for each updated design. It depends on 
 current temperature value of simulated annealing process.\cr
 This function is internally used by \code{\link{designGG}()}.
}
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
      either "SA" (simulated annealing) or "MH" 
  (Metropolis Hastings).
  }
  \item{temperature}{ 
      current temperature in simulated annealing process.
  }
}

\references{ 
  Y. Li, R. Breitling and R.C. Jansen "Generalizing genetical 
     genomics: the added value from environmental perturbation", Trends Genet
     (2008) 24:518-524. \cr
   E. Wit and J. McClure  "Statistics for Microarrays: Design, Analysis
      and Inference" (2004) Chichester: Wiley.
}
\author{  
 Yang Li <yang.li@rug.nl>, Gonzalo Vera <g.vera.rodriguez@rug.nl> \cr
 Rainer Breitling <r.breitling@rug.nl>, Ritsert Jansen <r.c.jansen@rug.nl> 
}
\keyword{method}
