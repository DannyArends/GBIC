\name{conditionLevel}
\alias{conditionLevel}
\title{ Levels of all environmental factors }
\description{
    Describe the levels of all environmental factors for each RIL in the
    experiment.\cr
    This is a subfunction needed for \code{designScore}, but is not directly used.
}
\usage{
    conditionLevel( array.allocation, condition.allocation, 
                    condition.combination,nEnvFactors )
}
\arguments{
  \item{array.allocation}{
      a matrix with nArray rows and nRIL columns.
      Elements of 1/0 indicates this RIL is/not selected 
      for this array.
  }
  \item{condition.allocation}{  
      a matrix with nCondition rows and nRIL columns.
      Elements of 1/0 indicates this RIL is/not selected 
      for this condition.
  }
  \item{condition.combination}{ 
      a matrix indicating all possible levels for environmental factors, 
      with dimension of nConditions by nEnvFactors.
  }
  \item{nEnvFactors}{
      number of environmental factors, integer.
      When \code{nEnvFactors} is bewteen 0 and 3, all main factor(s) and
      interacting facotr(s) will be included.
  }
}
\details{
    For one-color array, \code{array.allocation} is \code{NULL}. Then the \code{conditionLevel} is 
    decided by \code{condition.allocation}.
    For two-color array, \code{array.allocation} decides which RILs are selected and
    then the \code{condition.allocation} indicates which condition this RIL will be put in for the
    experiment. 
}
\value{
  A matrix with dimension of nRILs by nEnvFactors, each element indicates the 
  level of a certain environmental factor to which the RIL is exposed in the 
  experiment. 
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
  
\seealso{ \code{\link{designScore}}, \code{\link{conditionCombination}} }

\keyword{method}
