\name{updateDesign}
\alias{updateDesign}
\title{ Updates current design }
\description{
  Updates current experimental design (including \code{array.allocation} and \code{condition.allocation}).
}
\usage{
  updateDesign( array.allocation, condition.allocation, nRILs, 
                nSlides, nEnvFactors, nTuple, btwoColorArray )
}
\arguments{
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
  \item{nRILs}{ 
      number of RILs available for the experiment. 
  }
  \item{nSlides}{ 
      total number of slides available for experiment. 
  }
 \item{nEnvFactors}{
      number of environmental factors, integer.
      When \code{nEnvFactors} is between 0 and 3, all main factor(s) and
      interacting factor(s) will be included.
  }
   \item{nTuple}{
      average number of RILs to be assigned onto each condition. \cr
      \code{nTuple} should be a real number which is larger than 1.\cr
      If \code{nTuple} < 1, the algorithm will stop and show the message, \cr
      \code{warning: "The number of slides is too small to perform the experiment."}
  }
 \item{bTwoColorArray}{
      binary variable indicating experiment type: \cr
      \code{bTwoColorArray <- T} \#for dual channel experiment \cr
      \code{bTwoColorArray <- F} \#for single channel experiment 
  }
}
\details{
    This function calls two subfunctions: \code{conditionUpdate} and \code{arrayUpdate}.
}
\value{
  a list with two elements, \code{array.allocation} and \code{condition.allocation}.
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
\keyword{method}

 
