\name{pairLevel}
\alias{pairLevel}
\title{ Pair levels for paired RILs  }
\description{
    Pair levels for two RILs allocated into one slide(TwoColorArray)
    It is a subfunction needed for designScore function, but is not directly
    used.
}
\usage{
  pairLevel( xxx , rilNames )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xxx}{
     can be genotype.level, condition.level or interaction.level }
  \item{rilNames}{
     names for all RILs that have beedn selected for the experiment }

}
\details{
    This function is used only for two-color array.
}
\value{
  
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

\seealso{ See Also \code{\link{designScore}} }

\keyword{method}
