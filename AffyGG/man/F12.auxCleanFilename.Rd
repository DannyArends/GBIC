\name{auxCleanFilename}
\alias{auxCleanFilename}
\title{ Replace illegal characters in file names }
\description{
  This function replaces any illegal character i.e. ? or / )
  included in the input file name and returns a valid file name.
}
\usage{
auxCleanFilename( filename )
}
\arguments{
  \item{filename}{ a string with candidate file name. }
}
\details{
  In this implementation the following illegal characters:\cr
  \eqn{ ? [ ] / \\ = + < > : ; " , }\cr
  are replaced with the following ones:\cr
  \eqn{\# \{ \} \_ \_ - \_ \{ \} \_ \_ ' \_}    
}
\value{
  Example:
\tabular{lr}{
  Input: \tab traits.AFFX=YEL024w/RIP1\_at.csv  \cr
  Output: \tab traits.AFFX-YEL024w\_RIP1\_at.csv \cr
}     
}  
\author{ 
Gonzalo Vera <g.vera.rodriguez@rug.nl>
}
\note{
  This function doesn't handle relative or full path names. They should be 
  removed from the file name before using this function.\cr
  The illegal characters varies on each OS / file system. Here we try to 
  tackle the most common ones.\cr
  The user can easily modify the source code to map different characters. 
}
\seealso{ 
  \code{\link{splitProbeSetTraits}}
}
\examples{
## Not Run:

## For an example check source code of splitProbeSetTraits.R

## End(Not run)
}
