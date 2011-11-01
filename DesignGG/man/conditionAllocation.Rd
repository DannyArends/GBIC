\name{conditionAllocation}
\alias{conditionAllocation}
\title{ Allocate the selected RILs into different conditions
}
\description{
  This is a subfunction used by \code{initialDesign} but is not directly used.
  In  the experiment where samples are profiled in pairs, the samples are firstly
  selected and paired on each array and then the selected samples are randomly allocated
  into different conditions.
}
\usage{
  conditionAllocation( selectedRILs, genotype, nConditions, nSlides, nTuple )
}
\arguments{
  \item{selectedRILs}{ 
      the index of the selected RILs among all that are available for the experiment. 
  }
  \item{genotype}{
      genotype data: a nMarker-by-nRILs matrix with two allels being 0 and 1 (or A and B)
      or three allels being 0, 0.5 and 1 (or, A, H, and B), where 0.5 (or H) represents
      heterozygous allele.
  }
  \item{nConditions}{ 
      number of all possible combination of all environmental factors. It should be larger than 1. 
  }
  \item{nSlides}{
      total number of slides available for the experiment. It should be a non-zero integer.
  }
  \item{nTuple}{ 
      average number of RILs to be assigned onto each condition. \cr
      \code{nTuple} should be a real number which is larger than 1.\cr
      if \code{nTuple} < 1, the algorithm will stop and show a message as below, \cr
      \code{warning: "The number of slides is too small to perform the experiment."}
  }
}
\details{
    This function is only called by \code{initialDesign} function when \code{btwoColorArray}
     is \code{TRUE}.
}
\value{
  A matrix with nCondition rows and nRIL columns.
  Elements of 1/0 indicate that this RIL is/not selected for this condition.
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
\seealso{ \code{\link{initialDesign}}}

\keyword{method}