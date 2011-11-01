\name{updateDesign}
\alias{updateDesign}
\title{ Updates current design }
\description{
  Updates current design (including array.allocation and condition.allocation)
}
\usage{
  updateDesign( array.allocation, condition.allocation, nRILs, 
                nSlides, nEnvFactors, nTuple, btwoColorArray )
}
\arguments{
  \item{array.allocation}{ 
      matrix with nArray rows and nRIL columns.
      Elements of 1/0 indicates this RIL is/not selected 
      for this array. 
  }
  \item{condition.allocation}{ 
      matrix with nCondition rows and nRIL columns.
      elements of 1/0 indicates this RIL is/not selected 
      for this condition.  
  }
  \item{nRILs}{ 
      number of RILs available for the experiment 
  }
  \item{nSlide}{ 
      total number of slides available for experiment 
  }
  \item{nEnvFactors}{ 
      number of environmental factors, integer.
      When nEnvFactors is bewteen 0 and 3, all main factor and
      interacting facotrswill be included.
  }
  \item{nTuple}{ 
      average number of RILs put into each condition nTuple > 1, a real number.
      if nTuple < 1, stop,
      warning message:"The number slides is too less to perform the experiment."
  }
  \item{btwoColorArray}{ 
      binary variable indicating the array type in the experiment:
      bTwoColorArray = T, dual channel experiment
      bTwoColorArray = F, single channel experiment
  }
}
\details{
    This function calls two subfunctions: conditionUpdate and arrayUpdate
}
\value{
  a list with two elements, array.allocation and  condition.allocation
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
\keyword{method}

 
