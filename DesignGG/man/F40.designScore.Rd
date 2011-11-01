\name{designScore}
\alias{designScore}
\title{ calculate the A- or D- optimality score based on design}
\description{
    According to the current experiment design, the Fisher information matrix 
    is obtained and then either A- or D- optimality score is computed.
}
\usage{
  designScore( genotype, array.allocation, condition.allocation, 
               nEnvFactors, nLevels, Level, nConditions, weight=1, 
               optimality="A", bTwoColorArray, envFactorNames)
}
\arguments{
  \item{genotype}{
      genotype data: a matrix nMarker x nRILs with two allels being 0 and 1. }
  \item{array.allocation}{ 
      matrix iwth nArray rows and nRIL columns.
      elements of 1/0 indicates this RIL is/not selected 
      for this array.  
  }
  \item{condition.allocation}{ 
      matrix iwth nCondition rows and nRIL columns.
      elements of 1/0 indicates this RIL is/not selected 
      for this condition. 
  }
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
  \item{nConditions}{ 
      number of all possible combination of all environmental factors.  
  }
  \item{weight}{ 
       a vector with length of variablNumber which is calculated
      from function variablNumber, default = 1 which means the parameters to be 
      estimated are equally important during optimization.
  }
  \item{optimality}{ 
      type of optimality, i.e. "A" (A-optimality) or "D"
      (D-optimality). A-optimality minimizes $Trace((X'X)^{-1})$,
      which corresponds to minimum average variance of the
      parameter estimates. D-optimality minimizes $det(X'X)^{-1}$,
      which corresponds to minimum generalized variance of the
      parameter estimates.
  }
  \item{bTwoColorArray}{ 
      binary variable indicating experiment type:
      bTwoColorArray = T, dual channel experiment
      bTwoColorArray = F, single channel experiment
      }
  \item{envFactorNames}{
        a vector with names for all environmental factor(s). For example, for the
        experiment with two environmental factors of temperature and cell types
        envFactorNames <- c( "Temperature", "celltype" )
        default=NULL, then the output will use "var1" and "var2" to indicate the
        environmental factors.
  }
}
\details{
  Example parameter settings:
  Suppose to design an experiment with two environmental factors (F1, F2) and 
  there are two diffferent levels for each environment. The levels are 16
  and 24 for F1, and 5 anf 10 for F2. Thus the following command can be used:\cr
      > nEnvFactors 	<- 2       \cr
      > nLevels      	<- c ( 2, 2 ) \cr
      > levels 	      <- list ( c(16, 24), c(5, 10) )  \cr
  The length of parameter "weight" is dependent on the number of environmental
  factors: \cr
  When nEnvFactor = 0, \cr
    weight is 1 as there is only one parameter of interest (genotype).\cr
  When nEnvFactor = 1, \cr 
    weight = c( $w_{Q}$, $w_{F1}$, $w_{QF1}$ ) \cr
  When nEnvFactor = 2, \cr 
    weight = c( $w_{Q}$, $w_{F1}$, $w_{F2}$, $w_{QF1}$, $w_{QF2}$, $w_{F1F2}$, $w_{QF1F2}$) \cr
  When nEnvFactor = 3, \cr
    weight = c( $w_{Q}$, $w_{F1}$, $w_{F2}$, $w_{F2}$,
                $w_{QF1}$, $w_{QF2}$, $w_{QF3}$, $w_{F1F2}$, $w_{F1F3}$, $w_{F2F3}$,
                $w_{QF1F2}$, $w_{QF1F3}$, $w_{QF2F3}$, $w_{QF1F2F3}$ ) \cr    
  Here $w_{Q}$ represent the weight for genotype effect,  $w_{F1}$ represent the 
  weight for F1 effect and $w_{QF1}$ represent the weight for interaction between
  genotype and F1 effect, ect. 
}
\value{
  The score is defined as the "double" sum of the variances of summed over all 
  parameters and summed over all markers. 
}
\references{ 
   Y. Li, R. Breitling and R.C. Jansen "Generalizing genetical 
     genomics: the added value from environmental perturbation", Trends Genet
     (2008) 24:518-524. 
   E. Wit and J. McClure  "Statistics for Microarrays: Design, Analysis
      and Inference" (2004) Chichester: Wiley. 
  }
\author{  
 Yang Li <yang.li@rug.nl>, Gonzalo Vera <g.vera.rodriguez@rug.nl> \cr
 Rainer Breitling <r.breitling@rug.nl>, Ritsert Jansen <r.c.jansen@rug.nl> 
}

\seealso{ \code{\link{designGG}} }
\examples{}
\keyword{method}
