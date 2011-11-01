\name{createCytoFiles}
\alias{createCytoFiles}
\title{ Write Cytoscape visualization files }
\description{
Create visualization files for Cytoscape, the network visualization software.
Any correlation higher than \code{simThres} is written into Cytoscape network
files and can be loaded into Cytoscape for visualization.
}
\usage{
createCytoFiles(corrMatrix, filename, simThres = NULL, hideNodes = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{corrMatrix}{ matrix for pair-wise correlation.
  See \code{\link{corrSecondOrder}} example data. }
  \item{filename}{ network file name, without extension. Extensions .sif and .eda 
  will be added for network and edge attribute files, respectively. }
  \item{simThres}{ (optional) numeric similarity threshold if the values in 
  matrix are similarity values such as correlation coefficients. Default is NULL. }
  \item{hideNodes}{ (optional) logical value to hide nodes without any 
  significant links. Default is TRUE}
}
\value{
 A network file (\code{filename}.sif) and edge attribute file 
 (\code{filename}.eda) are generated.
}
\references{ 
  Fu J, Swertz MA, Keurentjes JJB, Jansen RC. MetaNetwork: a computational tool for the 
  genetic study of metabolism. Nature Protocols (2007). 
  
  \url{http://gbic.biol.rug.nl/supplementary/2007/MetaNetwork}
  \url{http://www.cytoscape.org/}
}
\author{
Jingyuan Fu <j.fu@rug.nl>, 
Morris Swertz <m.a.swertz@rug.nl>, 
Ritsert Jansen <r.c.jansen@rug.nl>
}
\seealso{ 
Use \code{\link{cor}}, \code{\link{qtlCorrZeroOrder}} and 
\code{\link{qtlCorrSecondOrder}} to calculate \code{corrMatrix} correlation matrix.
Use \code{\link{MetaNetwork}} for automated appliction of this function as 
part a genetic analysis protocol on metabolites.
}
\examples{
##NOTE: this method can be used on any correlation matrix.
#Here we use MetaNetwork methods.

## load the example data provided with this package                        
data(markers)   
data(genotypes)
data(traits)                                                       
                                             
##OR: load your own data                     
#markers        <- loadData("markers.csv")
#genotypes      <- loadData("genotypes.csv")
#traits         <- loadData("traits.csv")  
                                             
##calculate the two part qtl
qtlProfiles     <- qtlMapTwoPart(genotypes=genotypes, traits=traits[1:5,], spike=4)
  
##set the qtl threshold
qtlThres        <- 3.79

##OR: estimate the threshold yourself
#qtlThres       <- qtlThreshold(genotypes, traits, spike=4)  
  
##calculate zero order correlation
corrZeroOrder   <- qtlCorrZeroOrder(markers, qtlProfiles, qtlThres)

##calculate second order correlation
corrSecondOrder <- qtlCorrSecondOrder( corrZeroOrder ) 

##set the correlation threshold
corrThres       <- 0.14

##OR: estimate qtlCorrThreshold yourself
#corrThres      <- qtlCorThreshold(markers, genotypes, traits, spike=4, qtlThres=qtlThres)  

##create cytoscape files "mynetwork.sif" and "mynetwork.eda"
createCytoFiles(corrSecondOrder, "mynetwork", simThres = corrThres)

cat("cytofiles mynetwork.sif and mynetwork.eda created\n")
}