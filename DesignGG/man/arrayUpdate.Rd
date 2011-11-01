\name{arrayUpdate}
\alias{arrayUpdate}
\title{ Update array allocation }
\description{
    Update the allocation of samples on the arrays. 
    This is a subfunction needed for \code{updateDesign}, but is not directly used.
}
\usage{
  arrayUpdate(array.allocation, condition.allocation, nRILs, nSlides)
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
}
\details{
  This function is used only for designing a dual-channel experiment where samples
  need to be paired.
}
\value{
    A list with the following two elements: \cr
    \code{new.array.allocation}: an updated array allocation table \cr
    \code{new.condition.allocation}: an updated condition allocation table
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
\seealso{ 
  \code{\link{updateDesign}}
}
\keyword{method}
