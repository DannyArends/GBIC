\name{conditionAllocation}
\alias{conditionAllocation}
\title{ Allocate selected RILs into different conditions
}
\description{
  This is a subfunction used by initialDesign but is not directly used.
  In  the experiment where samples are profiled in pair, the samples were firstly
  selected and paired on each array and then the selected samples are allocated
  into different conditions.
}
\usage{
  conditionAllocation( selectedRILs, genotype, nConditions, nSlides, nTuple )
}
\arguments{
  \item{selectedRILs}{ 
      selected RILs from the array.allocation 
  }
  \item{genotype}{
      genotype data: a matrix nMarker x nRILs with two allels being 0 and 1.
  }
  \item{nConditions}{ 
      number of all possible combination of all environmental factors 
  }
  \item{nSlides}{
      total number of slides available for experiment
  }
  \item{nTuple}{ 
      average number of RILs put into each condition \cr
      nTuple > 1, a real number.\cr
      if nTuple < 1, stop, \cr
        \code{warning: "The number slides is too less to perform the experiment."}
  }
}
\details{
    This function is only called by initialDesign function when btwoColorArray
     is \code{TRUE}.
}
\value{
  A matrix with nCondition rows and nRIL columns.
  Elements of 1/0 indicates this RIL is/not selected for this condition.
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
\seealso{ \code{\link{initialDesign}}}

\keyword{method}
