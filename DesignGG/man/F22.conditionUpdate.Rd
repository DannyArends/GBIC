\name{conditionUpdate}
\alias{conditionUpdate}
\title{ Update condition allocation }
\description{
    This is a subfunctions needed for updateDesign, but is not directly used.
}
\usage{
    conditionUpdate( condition.allocation, nTuple )
}
\arguments{
  \item{condition.allocation}{ 
      matrix iwth nCondition rows and nRIL columns.
      elements of 1/0 indicates this RIL is/not selected 
      for this condition. 
  }
  \item{nTuple}{ 
      average number of RILs put into each condition 
      nTuple > 1, a real number. \cr
      if nTuple < 1, stop, \cr
        warning: "The number slides is too less to perform the experiment."
  }
}
\details{
    This function will be used both in single and dual channel experiment 
    design.
}
\value{
  An updated condition.allocation table.
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

\seealso{ \code{\link{arrayUpdate}}, \code{\link{designGG}}  }
\keyword{method}
