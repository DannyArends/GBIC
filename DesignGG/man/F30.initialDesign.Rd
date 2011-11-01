\name{initialDesign}
\alias{initialDesign}
\title{ Initilize a design 
 }
\description{
 allocate RILs into different conditional and pair RILs into sldies

}
\usage{
  initialDesign( genotype, nRILs, nSlides, nConditions, nTuple, 
                 bTwoColorArray )
}

\arguments{
  \item{genotype}{ 
      genotype data: a matrix nMarker x nRILs with two allels 
      being 0 and 1 
  }
  \item{nRILs}{ 
      total number of RILs available for experiment 
  }
  \item{nSlides}{ 
      total number of slides available for experiment 
  }
  \item{nConditions}{ 
      number of all possible combination of all environmental factors  
  }
  \item{nTuple}{  
      average number of RILs put into each condition 
      nTuple > 1, a real number. \cr
      if nTuple < 1, stop,
      warning: "The number slides is too less to perform the experiment."
  }
  \item{bTwoColorArray}{ 
      binary variable indicating the array type in the experiment: \cr
      bTwoColorArray = T, dual channel experiment \cr
      bTwoColorArray = F, single channel experiment
  }
}
\details{
    For two-color array, randomly choose a RIL and pair it with its 
    mostly genetically different RIL on one array. \cr
    For one-color array experimen, array.allocation is NULL as there is no 
    need to pair samples.
}
\value{
  a list with 2 matrices:  \cr
  condition.allocation: allocate RILs into different conditional
    (nrow = nCondition, ncol= nRILs)      \cr
  array.allocation: pair RILs into sldies
    (nrow = nSlide, ncol= nRILs) 
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
\note{ This function calls conditionAllocation to allocate selected RILs 
            into different conditions.
}
\seealso{ \code{\link{designGG}} }
\examples{
	> genotype     	<- read.csv ("genotype.csv")
	> nEnvFactors 	<- 2  
  > nLevels      	<- c ( 2, 2 )
  > levels 	      <- list ( c(16, 24), c(5, 10) )
  > nSlides  	    <- 100
  > nTuple	      <- 25  
  > bTwoColorArray<- TRUE 
  >initialDesign( genotype, nRILs, nSlides, nConditions, nTuple, 
                 bTwoColorArray )          
 }
\keyword{method}
