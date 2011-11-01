\name{variableNames}
\alias{variableNames}
\title{ Generate variable names for all factors }
\description{
  Generate variable names for genetic, environmental factors and interacting
        terms
}
\usage{
variableNames(nEnvFactors)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nEnvFactors}{ 
      number of environmental factors, integer.
      When nEnvFactors is bewteen 0 and 3, all main factor and
      interacting facotrswill be included.
  }
}
\details{
  generates names for variables, length = variableNumber+1
}
\value{
   When nEnvFactor = 0, variableNumber = 1  \cr
   variableNames = one genetic factor "Q".  \cr
   When nEnvFactors = 1, variableNumber = 3  \cr
   variableNames = one genetic factor "Q",    \cr
                   one envrinmental factor "F", \cr
                   and one interacting factor "QxF". \cr   
  When nEnvFactor = 2, \cr
  variableNames = one genetic factor "Q", \cr
                  two envrinmental factors "F1" and "F2", \cr
                  three two-way interacting factors "QF1", "QF2", "F1F2", \cr
                  and one three way interacting factors "QxF1xF2".  \cr
  When nEnvFactor = 3, \cr
  variableNames = one genetic factor "Q", \cr
                  three envrinmental factors "F1", "F2" and "F3", \cr
                  six two-way interacting factors "QF1", "QF2", "QF3", "F1F2",
                   "F2F3" and "F1F3",\cr
                  four three-way interacting factors "QxF1xF2", "QxF1xF3",
                   "QxF2xF3", "F1xF2xF3" \cr
                  and one four-way interacting factors "QxF1xF2xF3". \cr

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

\seealso{ See Also \code{\link{variableNumber} }  }
\keyword{method}
