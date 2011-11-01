\name{qtlCorrSecondOrder}
\alias{qtlCorrSecondOrder}
\title{ Calculate the second-order correlation  }
\description{
  Calculate the second-order partial correlation between any pair of QTL profiles to
  decompose the confounded effect of multiple QTLs.
}
\usage{
qtlCorrSecondOrder(corrZeroOrder, topCorNo=20, filename = NULL)
}
\arguments{
  \item{corrZeroOrder}{ a zero order correlation matrix that has been calculated before.
  See \code{\link{corrZeroOrder}} example data.}
  \item{topCorNo}{ (optional) to reduce computing time, limit the pairwise 
  conditioning to only the \code{topCorNo} top correlated compounds. 
  Default is 20. To condition on all pairs, set topCorNo to nrow(corrZeroOrder).}
  \item{filename}{ (optional) path of the file where the correlations are to be stored. Default NULL. }
}
\details{
First-order correlation between \eqn{x} and \eqn{y} conditional on a single 
variable \eqn{z} are calculated by

        \deqn{r_{xy \arrowvert z} = \frac {r_{xy} - r_{xz}r_{yz}}
        {\sqrt{(1-r_{xz}^2)(1-r_{yz}^2)}}} 


where \eqn{r_{xy}}, \eqn{r_{xz}} and \eqn{r_{yz}} are zero-order correlation coefficients 
between traits \eqn{x} and \eqn{y}, \eqn{x} and \eqn{z}, 
and \eqn{y} and \eqn{z}, respectively (provided as \code{corrZeroOrder}. 

Then second-order partial correlation between \eqn{x} and \eqn{y}, conditional on a pair 
of variables \eqn{z} and \eqn{k} is a function of the first-order coefficients.

  \deqn{r_{xy\arrowvert zk} = \frac{r_{xy\arrowvert z} - r_{xk\arrowvert z}r_{yk\arrowvert z}}
        {\sqrt{(1-r_{xk\arrowvert z}^2)(1-r_{yk\arrowvert z}^2)}}}
                
For each pair \eqn{x} and \eqn{y}, the second-order partial correlations are calculated conditional
on each pair \eqn{z} and \eqn{k} and the minimal value is stored.

In order to save computing time, conditional variables z and k can 
be chosen from \code{topCorNo} top correlated variable x and y (e.g. \code{topCorNo} = 20).
}
\value{
    Returns a matrix of second-order partial correlations.
}
\source{  
 Keurentjes JJB, FU J, de vos CHR, Lommen A, Hall RD, Bino RJ, van der Plas LHW, Jansen RC, 
 Vreugdenhil D, and Koornneef M. The genetics of plant metabolism. Nature Genetics (2006) 7: 842-849.   
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
Use \code{\link{qtlCorrZeroOrder}} to calculate \code{corrZeroOrder} zero order correlation.
Use \code{\link{MetaNetwork}} for automated appliction of this function as part 
a genetic analysis protocol on metabolites.
}
\examples{
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

##calculate the second order correlation
corrSecondOrder <- qtlCorrSecondOrder(corrZeroOrder=corrZeroOrder)

##view the correlations
corrSecondOrder[1:5,1:5]
}