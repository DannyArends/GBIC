\name{pairLevel}
\alias{pairLevel}
\title{ Pair levels for paired RILs  }
\description{
    Pair levels for two RILs allocated into one slide (\code{bTwoColorArray=T}). It 
    is a subfunction needed for \code{designScore} function, but is not directly
    used.
}
\usage{
  pairLevel( xxx, rilNames )
}

\arguments{
  \item{xxx}{
     can be \code{genotype.level}, \code{condition.level} or \code{interaction.level} }
  \item{rilNames}{
     names for all RILs that have been selected for the experiment }

}
\details{
    This function is used only for two-color array.
}
\value{TODO}
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

\seealso{ See Also \code{\link{designScore}} }

\keyword{method}
