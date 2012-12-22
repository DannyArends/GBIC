\name{probeElimination}
\alias{probeElimination} 
\title{Probe elimination procedure}

\description{
  Performs the probe backwards elimination procedure per
  probe set. Probes are backwards eliminated to reduce the probe \eqn{\times} QTL
  interaction effect.
}

\usage{
probeElimination( probesetName, markerName, genotypes, traits, 
                  batch, filename=NULL, appendRecord=F )
}

\arguments{

  \item{probesetName}{ name of the probe set. }    
  \item{markerName}{ name of the marker. }
  \item{genotypes}{ matrix of genotypes for each marker (rownames) and individual 
  (columnnames), as numeric values 1, 2 or NA when missing.\cr
  See \code{\link{genotypes}} example data.}
  \item{traits}{ matrix of phenotypes for each trait (rownames) and individual 
  (columnnames), as numeric or NA when missing.\cr
  See \code{\link{traits}} and \code{\link{traits2}} example data. }
  \item{batch}{ vector with batch numbers in which each individual (e.g. BXD)
  was processed.} 
  \item{filename}{ (optional) path of the file where the data will be loaded 
  from or save to.\cr
  In this function, output file where results are stored (*.csv = comma
  separated values). Default is NULL. }  
  \item{appendRecord}{ (optional) boolean value to indicate whether results are appended 
  to an existing file (\code{T} or \code{F}). Default is \code{F}. }
}

\details{
  Performs the probe backwards elimination procedure per probe set for the cis 
  acting marker. Probes are backwards eliminated to reduce the probe \eqn{\times}{x}
  QTL interaction effect. Probes are eliminated until 2 probes are left. Using the wished 
  threshold, the user can find in the results which probes were eliminated and 
  in which order. 
  The results also show whether the first allele had higher expression than the second
  allele (i.e. \code{T}) in the eliminated probes and in the non-eliminated probes.
}

\value{
  Returns a \code{data.frame} with the following columns:

 \item{dropped}{ the number of the probe being eliminated. } 
 \item{pmarker}{ QTL effect (P value) of the remaining probes. }   
 \item{pinteraction}{ QTL \eqn{\times} probe effect (P value) of the remaining probes. }  
 \item{pmarkerdropped}{ QTL effect (P value) of the dropped probes.}  
 \item{pintdropped}{ QTL \eqn{\times} probe effect (P value) of the dropped probes. }  
 \item{alarm}{ \code{TRUE} if QTL effect dropped probes \eqn{>} 
                              QTL effect remaining probes.}  
 \item{allele1bigger}{ \code{TRUE} if first allele has higher signals than 
                        second allele in remaining probes. }  
 \item{allele1biggerdropped}{ \code{TRUE} if first allele has higher signals than 
                        second allele in dropped probes.} 
 Example table:
\tabular{lrrrrr}{
  \tab dropped \tab pmarker \tab pinteraction \tab pmarkerdropped    \tab pintdropped \cr
1 \tab 0  \tab 0.000000e+00 \tab 3.743150e-153   \tab 0.000000e+00   \tab 0.000000e+00 \cr 
2 \tab 15 \tab 2.078871e-09 \tab 7.602871e-121  \tab -1.000000e+00  \tab -1.000000e+00 \cr
3 \tab 16 \tab 1.723181e-02  \tab 1.146263e-03   \tab 1.654960e-15   \tab 5.356472e-03 \cr
4 \tab 3  \tab 6.480820e-03  \tab 9.393928e-03   \tab 9.107300e-15   \tab 2.407649e-28 \cr
5 \tab 9  \tab 6.679409e-03  \tab 5.115019e-02   \tab 1.576042e-14   \tab 2.037995e-38 \cr
}
\tabular{lrrr}{
 \tab alarm \tab allele1bigger \tab allele1biggerdropped \cr
 \tab FALSE         \tab FALSE                \tab FALSE \cr 
 \tab TRUE          \tab TRUE                \tab FALSE \cr
 \tab TRUE          \tab TRUE                 \tab TRUE \cr
 \tab TRUE          \tab TRUE                 \tab TRUE \cr
 \tab TRUE          \tab TRUE                 \tab TRUE \cr
}
}

\references{
Alberts R, Terpstra P, Li Y, Breitling R, Nap JP, Jansen RC.
Sequence polymorphisms cause many false cis eQTLs. 
\eqn{Submitted}
}

\author{ 
Rudi Alberts <r.alberts@rug.nl>, 
Ritsert Jansen <r.c.jansen@rug.nl>
}

\seealso{
   \code{\link{rma.preprocessing}},  \code{\link{qtlMap.xProbeSet}}
}

\examples{
# We load and define the data needed
data( traits2 ) 
data( genotypes )
data( batch )

# We call the probe elimination function
pe <- probeElimination( "160371_at", "D7Mit100", genotypes, traits2, batch )
}


