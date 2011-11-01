\name{conditionCombination}
\alias{conditionCombination}
\title{ Generate a matrix indicating all possible levels for environmental 
    factors
 }
\description{
    Generate a matrix indicating all possible levels for environmental 
    factors with dimension nConditions x nEnvFactors.
    This is a subfunctions needed for designScore, but is not directly used.
}
\usage{
  conditionCombination( nEnvFactors, nLevels, Level, envFactorNames )
}
\arguments{
  \item{nEnvFactors}{ 
      number of environmental factors, integer.
      When nEnvFactors is bewteen 0 and 3, all main factor and
      interacting facotrswill be included. 
  }
  \item{nLevels}{ 
      number of levels for each factor, a vector with each 
      component being integer. The length of it should equal 
      nEnvFactors.
  }
  \item{Level}{ 
      a list which specifies the levels for each factor in the
      experiment. There are in total nEnvFactors elements in the list and each 
      element correpsond to certain envrironmental factor. The emlemet is a 
      vector describing all levels of the environmental factor.     
      default setting for the level of each factor is 1,2,....nlevels[i]. (Here
      nLevels[i] is the ith element of nLevels, which tells the total number of 
      levels for ith environmental facotor) 
  }
  \item{envFactorNames}{
        a vector with names for all environmental factor(s). For example, for the
        experiment with two environmental factors of temperature and cell types
        envFactorNames <- c( "Temperature", "celltype" ) \cr
        default=NULL, then the output will use "var1" and "var2" to indicate the
        environmental factors.
  }
}
\details{
    Currently this function works only for nEnvFactors is between 0 and 3.
}

\value{
  a matrix with dimension of nConditions x nEnvFactors.
  Each element in the matrix indicates the level of corresponding environmental
   factor.

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

\seealso{ \code{\link{designScore}}}
\keyword{method}
