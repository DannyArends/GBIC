\name{variableNames}
\alias{variableNames}
\title{ Generate variable names for all factors }
\description{
  Generate variable names for genetic, environmental factors and interacting
        terms.
}
\usage{
variableNames(nEnvFactors)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{nEnvFactors}{
      number of environmental factors, integer.
      When \code{nEnvFactors} is bewteen 0 and 3, all main factor(s) and
      interacting facotr(s) will be included.
  }
}
\details{
  generates names for variables, a vector with the length of  (variableNumber+1).
}
\value{
   When \code{nEnvFactors} = 0, \code{variableNumber} = 1, and
   \code{variableNames} is one genetic factor "Q".  \cr
   When \code{nEnvFactors} = 1, \code{variableNumber} = 3, and
   \code{variableNames} are one genetic factor "Q", one environmental factor "F",
    and one interacting factor "QxF". \cr   
  When \code{nEnvFactors} = 2, \code{variableNumber} = 7, and \code{variableNames} are one genetic factor "Q",
                  two environmental factors "F1" and "F2", 
                  three two-way interacting factors "QF1", "QF2", "F1F2", 
                  and one three way interacting factors "QxF1xF2".  \cr
  When \code{nEnvFactors} = 3,  \code{variableNumber} = 15, and
  \code{variableNames} are one genetic factor "Q",
                  three environmental factors "F1", "F2" and "F3",
                  six two-way interacting factors "QF1", "QF2", "QF3", "F1F2",
                   "F2F3" and "F1F3",
                  four three-way interacting factors "QxF1xF2", "QxF1xF3",
                   "QxF2xF3", "F1xF2xF3" 
                  and one four-way interacting factors "QxF1xF2xF3". 
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

\seealso{ \code{\link{variableNumber} }  }
\keyword{method}
