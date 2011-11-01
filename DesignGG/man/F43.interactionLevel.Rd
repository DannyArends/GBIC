\name{interactionLevel}
\alias{interactionLevel}
\title{ Generate levels for all interacting factors }
\description{
    generate levels for all interacting factors for all RILs
    This is a subfunctions needed for designScore, but is not directly used.
}
\usage{
  interactionLevel( genotype.level, condition.level, markerIndex, 
                    nEnvFactors )
}
\arguments{
  \item{genotype.level}{ 
      levels of genetic factor for each RIL in the experiment  
  }
  \item{condition.level}{ 
      levels of all environmental factors for each RIL in the experiment 
  }
  \item{markerIndex}{ 
      indicate which genome position that level of genetic factor corresponds to 
  }
  \item{nEnvFactors}{ 
      number of environmental factors, integer.
      When nEnvFactors is bewteen 0 and 3, all main factor and
      interacting facotrswill be included.
  }
}
\details{
    markerIndex indicates the genome position that genotype.level correspond
    to. Experiment design is defined to be optimal over all markers if 
    the sum of scores, e.g. A-optimality criterion over all markers is 
    minimized.
}
\value{
  a matrix with nRILs rows. The number columns depends on the nEnvFactors.
  For example: \cr
    If nEnvFactors = 1, there is only one interaction term. \cr
    If nEnvFactors = 2, there is 3 pair-wise two-way interaction and 1 
        three-way interaction term.
}
\references{ 
 Y. Li, R. Breitling and R.C. Jansen "Generalizing genetical 
     genomics: the added value from environmental perturbation", Trends Genet
     (2008) 24:518-524. 
}
\author{  
 Yang Li <yang.li@rug.nl>, Gonzalo Vera <g.vera.rodriguez@rug.nl> \cr
 Rainer Breitling <r.breitling@rug.nl>, Ritsert Jansen <r.c.jansen@rug.nl> 
}

\seealso{ \code{\link{designScore}} ,\code{\link{conditionLevel}}  }

\keyword{method}
