\name{designScore}
\alias{designScore}
\title{ Calculate the A- or D- optimality score based on current experimental design}
\description{
    According to the current experimental design, the Fisher information matrix 
    is obtained and then either the A- or D- optimality score is computed.
}
\usage{
  designScore( genotype, array.allocation, condition.allocation, 
               nEnvFactors, nLevels, Level, nConditions, weight=1, 
               optimality="A", bTwoColorArray, envFactorNames)
}
\arguments{
  \item{genotype}{
      genotype data: a nMarker-by-nRILs matrix with two allels being 0 and 1 (or A and B)
      or three allels being 0, 0.5 and 1 (or, A, H, and B), where 0.5 (or H) represents
      heterozygous allele.
  }
 \item{array.allocation}{ 
      matrix with nArray rows and nRIL columns.
      Elements of 1/0 indicate this RIL is/not selected 
      for this array. 
  }
   \item{condition.allocation}{ 
      matrix with nCondition rows and nRIL columns.
      Elements of 1/0 indicate this RIL is/not selected 
      for this condition. 
  }
  \item{nEnvFactors}{
      number of environmental factors, integer.
      When \code{nEnvFactors} is between 0 and 3, all main factor(s) and
      interacting facotr(s) will be included.
  }
 \item{nLevels}{
      number of levels for each factor, a vector with each
      component being an integer. The length of it should equal
      \code{nEnvFactors}.
  }
\item{Level}{
      a list which specifies the levels for each factor in the
      experiment. There are in total \code{nEnvFactors} elements in the list and each
      element correpsond to certain envrironmental factor. The emlemet is a
      vector describing all levels of the environmental factor.
      default setting for the level of each factor is 1, 2, ... nLevels[i]. (Here
      nLevels[i] is the \emph{i}th element of nLevels, which gives the total number of
      levels for \emph{i} environmental factor).
  }
  \item{nConditions}{ 
      number of all possible combination of all environmental factors.  
  }
  \item{weight}{
      a vector with length of variableNumber which is calculated
      from function \code{variableNumber}. Default = 1 (which means the parameters to be 
      estimated are equally important during optimization.)
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
      binary variable indicating experiment type: \cr
      \code{bTwoColorArray <- T} \#for dual channel experiment \cr
      \code{bTwoColorArray <- F} \#for single channel experiment 
  }
   \item{envFactorNames}{
        a vector with names for all environmental factor(s). For example, for the
        experiment with two environmental factors of temperature and drug treatment:
        \code{envFactorNames <- c( "Temperature", "Dosage" ) }\cr
        Default = \code{NULL}, then the output will use "F1" and "F2" to indicate the
        environmental factors.
  }
}
\details{
  Example parameter settings:\cr
  Suppose to design an experiment with two environmental factors (F1, F2) and 
  there are two diffferent levels for each environment. The levels are 16
  and 24 for F1, and 5 anf 10 for F2. Thus the following command can be used:\cr
  \code{nEnvFactors <- 2}\cr
  \code{nLevels <- c ( 2, 2 )}\cr
  \code{levels <- list ( c(16, 24), c(5, 10) )}\cr
  The length of parameter \code{weight} is dependent on the number of environmental
  factors: \cr
  When \code{nEnvFactor} = 0, \cr
    \code{weight} is 1 as there is only one parameter of interest (genotype).\cr
  When \code{nEnvFactor} = 1, \cr 
    \code{weight} = c( $w_{Q}$, $w_{F1}$, $w_{QF1}$ ) \cr
  When \code{nEnvFactor} = 2, \cr 
    \code{weight} = c( $w_{Q}$, $w_{F1}$, $w_{F2}$, $w_{QF1}$, $w_{QF2}$, $w_{F1F2}$, $w_{QF1F2}$) \cr
  When \code{nEnvFactor} = 3, \cr
    \code{weight} = c( $w_{Q}$, $w_{F1}$, $w_{F2}$, $w_{F2}$,
                $w_{QF1}$, $w_{QF2}$, $w_{QF3}$, $w_{F1F2}$, $w_{F1F3}$, $w_{F2F3}$,
                $w_{QF1F2}$, $w_{QF1F3}$, $w_{QF2F3}$, $w_{QF1F2F3}$ ) \cr    
  Here $w_{Q}$ represents the weight for genotype effect,  $w_{F1}$ represent the 
  weight for F1 effect and $w_{QF1}$ represent the weight for interaction between
  genotype and F1 effect, etc. 
}
\value{
  The score is defined as the "double" sum of the variances, summed over all 
  parameters and over all markers. 
}
\references{ 
  Y. Li, M. Swertz, G. Vera, J. Fu, R. Breitling, and R.C. Jansen. designGG: 
  An R-package and Web tool for the optimal design of genetical genomics 
  experiments. (submitted) \cr
  http://gbic.biol.rug.nl/designGG \cr
  Y. Li, R. Breitling and R.C. Jansen. Generalizing genetical 
     genomics: the added value from environmental perturbation, Trends Genet
     (2008) 24:518-524. \cr
  E. Wit and J. McClure. Statistics for Microarrays: Design, Analysis
      and Inference. (2004) Chichester: Wiley.    
  }
\author{ 
Yang Li <yang.li@rug.nl>, Gonzalo Vera <gonzalo.vera.rodriguez@gmail.com> \cr
  Rainer Breitling <r.breitling@rug.nl>, Ritsert Jansen <r.c.jansen@rug.nl> 
}

\seealso{ \code{\link{designGG}} }
\keyword{method}
