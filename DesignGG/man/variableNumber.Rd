\name{variableNumber}
\alias{variableNumber}
\title{ Compute the number of variables in the experiment }
\description{
    \code{nEnvFactors} = 0, only genetic factor is considered. \cr
    \code{nEnvFactors} > 1, genetic and environmental facotrs, and all possible
    interacting factors are considered. \cr     
}
\usage{
  variableNumber( nEnvFactors )
}
\arguments{
  \item{nEnvFactors}{
      number of environmental factors, an integer.\cr
      When nEnvFactors is between 0 and 3, all main factors and
      interacting factors will be included.
  }
}

\value{
  \code{nEnvFactors} = 1, \code{variableNumber} = 3 (one genetic factor Q, one environmental 
                                    factor F, and one interacting factor QxF)\cr  
  \code{nEnvFactors} = 2, \code{variableNumber} = 7   \cr
  \code{nEnvFactors} = 3, \code{variableNumber} = 15   \cr
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

\seealso{ \code{\link{variableNames}}}
\examples{
#TODO
}
\keyword{method}
