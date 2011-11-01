\name{interactionLevel}
\alias{interactionLevel}
\title{ Generate levels for all interacting factors }
\description{
    Generate levels for all interacting factors for all RILs.
    This is a subfunction needed for \code{designScore}, but is not directly used.
}
\usage{
  interactionLevel( genotype.level, condition.level, markerIndex, 
                    nEnvFactors )
}
\arguments{
  \item{genotype.level}{ 
      levels of genetic factor for each RIL in the experiment.  
  }
  \item{condition.level}{ 
      levels of all environmental factors for each RIL in the experiment. 
  }
  \item{markerIndex}{ 
      indicate which genome position that level of genetic factor corresponds to. 
  }
  \item{nEnvFactors}{
      number of environmental factors, integer.
      When \code{nEnvFactors} is bewteen 0 and 3, all main factor(s) and
      interacting facotr(s) will be included.
  }
}
\details{
    \code{markerIndex} indicates the genome position that \code{genotype.level} corresponds
    to.\cr 
    An experiment design is defined to be optimal over all markers if 
    the sum of scores, e.g. A-optimality criterion over all markers is 
    minimized.
}
\value{
  a matrix with nRILs rows. The number columns depends on \code{nEnvFactors}.
  For example: \cr
    If \code{nEnvFactors} = 1, there is only one interaction term. \cr
    If \code{nEnvFactors} = 2, there are three pair-wise two-way interaction terms and one 
        three-way interaction term.
}
\references{ 
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
  Rainer Breitling <r.breitling@rug.nl>, Ritsert Jansen <r.c.jansen@rug.nl> 
}

\seealso{ \code{\link{designScore}}, \code{\link{conditionLevel}}  }

\keyword{method}
