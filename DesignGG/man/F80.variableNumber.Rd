\name{variableNumber}
\alias{variableNumber}
\title{ Compute the number of variables in the experiment }
\description{
    nEnvFactors=0, only genetic factor is considered. \cr
    nEnvFactors>1, genetic factor, envrinmental facotrs, and all possible
    interacting factors are considered. \cr     
}
\usage{
  variableNumber( nEnvFactors )
}
\arguments{
  \item{nEnvFactors}{
      number of environmental factors, integer.\cr
      When nEnvFactors is bewteen 0 and 3, all main factor and
      interacting facotrswill be included.
  }
}

\value{
  nEnvFactors = 1, variableNumber = 3 (one genetic factor Q, one envrinmental 
                                    factor F, and one interacting factor QxF)\cr  
  nEnvFactors = 2, variableNumber = 7   \cr
  nEnvFactors = 3, variableNumber = 15   \cr
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

\seealso{ See Also \code{\link{variableNames}}}
\examples{}
\keyword{method}
