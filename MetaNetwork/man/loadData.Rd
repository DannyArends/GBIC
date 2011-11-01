\name{loadData}
\alias{loadData}
\title{ Load data from csv files }
\description{
  Load a comma separated values (csv) file and convert the contents of the file 
  into a matrix object, taking the 1st column as the rownames and the 1st
  row as columnnames. See also \code{\link{read.csv}}.
}
\usage{
  mydata <- loadData(filename)
}
\arguments{
  \item{filename}{ path of the file where the data will be loaded from. }
}
\value{
  Returns a matrix containing the data from the indicated file.
}
\references{ 
  Fu J, Swertz MA, Keurentjes JJB, Jansen RC. MetaNetwork: a computational tool for the 
  genetic study of metabolism. Nature Protocols (2007). 
  
  \url{http://gbic.biol.rug.nl/supplementary/2007/MetaNetwork} 
}
\author{
Jingyuan Fu <j.fu@rug.nl>, 
Morris Swertz <m.a.swertz@rug.nl>, 
Ritsert Jansen <r.c.jansen@rug.nl>
}
\seealso{ 
Use \code{\link{genotypes}}, \code{\link{traits}} and \code{\link{markers}} 
as example data sets.
Use \code{\link{MetaNetwork}} for automated analysis of loaded data sets as 
part of a genetic analysis protocol on metabolites.
}
\examples{
  \dontrun{
  ##load genotypes
  genotypes <- loadData("genotypes.csv") 
  
  ##load traits  
  traits    <- loadData("traits.csv")
  
  ##load markers   
  markers   <- loadData("markers.csv")
  }
}