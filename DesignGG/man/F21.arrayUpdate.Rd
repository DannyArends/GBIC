\name{arrayUpdate}
\alias{arrayUpdate}
\title{ Update array allocation }
\description{
    DESCRIPTION MISSING. \cr
    This function is internally used by \code{\link{updateDesign}()}.
}
\usage{
  arrayUpdate(array.allocation, condition.allocation, nRILs, nSlides)
}
\arguments{
  \item{array.allocation}{ 
      matrix with \code{nArray} rows and \code{nRIL} columns.
      elements of 1/0 indicates this RIL is/not selected 
      for this array. 
  }
  \item{condition.allocation}{ 
      matrix with \code{nCondition} rows and \code{nRIL} columns.
      elements of 1/0 indicates this RIL is/not selected 
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
  This function is used only for designing a dual channel experiment where samples
  need to be paired.
}
\value{
    A list with two elements: \cr
    \item{new.array.allocation}{ an updated array allocation table. }
    \item{new.condition.allocation}{ an updated condition allocation table.}
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
\seealso{ 
  \code{\link{updateDesign}}
}
\keyword{method}
