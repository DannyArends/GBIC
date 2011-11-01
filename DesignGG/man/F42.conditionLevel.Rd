\name{conditionLevel}
\alias{conditionLevel}
\title{ levels of all environmental factors }
\description{
    Describe the levels of all environmental factors for each RILs in the
    experiment.\cr
    This is a subfunctions needed for designScore, but is not directly used.
}
\usage{
    conditionLevel( array.allocation, condition.allocation, 
                    condition.combination,nEnvFactors )
}
\arguments{
  \item{array.allocation}{
      matrix with nArray rows and nRIL columns.
      Elements of 1/0 indicates this RIL is/not selected 
      for this array.
  }
  \item{condition.allocation}{  
      matrix with nCondition rows and nRIL columns.
      Elements of 1/0 indicates this RIL is/not selected 
      for this condition.
  }
  \item{condition.combination}{ 
      a matrix indicating all possible levels for environmental factors, 
      with dimension of nConditions by nEnvFactors.
  }
  \item{nEnvFactors}{
      number of environmental factors, integer.
      When nEnvFactors is bewteen 0 and 3, all main factor and
      interacting facotrswill be included.
  }
}
\details{
    For one-color array, array.allocation is NULL. Then the conditionLevel is 
    decided by condition.allocation.
    For two-color array, array.allocation decides which RILs are selected and
    then the condition.allocation gives which condition this RIL is put in the
    experiment. 
}
\value{
  A matrix with dimension of nRILs x nEnvFactors, each element indicate the 
  level of certain environmental factor at which the RIL is exposed in the 
  experiment. 
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
  
\seealso{ \code{\link{designScore}}, \code{\link{conditionCombination}} }

\keyword{method}
