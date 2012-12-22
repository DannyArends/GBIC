\name{rma.preprocessing}
\alias{rma.preprocessing}
\title{Preprocessing Affymetrix .CEL files using RMA}

\description{
  Applies the RMA preprocessing method to Affymetrix 
  .CEL files. Output is obtained on probe level (without summarizing).
}

\usage{
rma.preprocessing( directory, celfiles, filename, 
                   cdffile, cdfpackage )
}

\arguments{
  \item{directory}{ path of the directory where files are located.\cr
  In this function, input directory where .CEL files are located. }

  \item{celfiles}{ name of the .CEL files that should be processed. }  

  \item{filename}{ (optional) path of the file where the data will be loaded from or saved to.\cr
  In this function, output file where results are stored (*.csv = comma
  separated values). Default is NULL.}
  
  \item{cdffile}{ (optional) filename argument with a custom Affymetrix chip 
                   description file (CDF). Default is NULL.}
  
  \item{cdfpackage}{ (optional) filename argument with a package containing the CDF file name. 
                   Default is NULL.}
}

\details{
  \code{\link{rma.preprocessing}} does RMA (Robust Multichip Analysis) background 
  correction and quantile normalization. Probe level results are returned 
  and optionally stored. 
}

\value{
  Returns a \code{data.frame} (and creates a .csv file) 
  where on each line the probe signal values for all .CEL files supplied are given. \cr
  Example with 3 .CEL files:
  \preformatted{
> celfiles <- c("bxd1a.cel","bxd2a.cel","bxd5a.cel")    
> probesignals <- rma.preprocessing( directory, celfiles, 
                                     "probesignals.csv" )
> dim(probesignals)
[1] 197993      4
> head(probesignals)}\cr
\tabular{lrrrr}{
  \tab probeset \tab bxd1a.cel \tab bxd2a.cel \tab bxd5a.cel \cr
1 \tab 100001\_at     \tab 4.935     \tab 4.840     \tab 4.950 \cr
2 \tab 100001\_at    \tab 10.514    \tab 10.610    \tab 10.385 \cr
3 \tab 100001\_at     \tab 9.396     \tab 9.849     \tab 9.651 \cr
4 \tab 100001\_at    \tab 10.682    \tab 11.427    \tab 10.482 \cr
5 \tab 100001\_at    \tab 11.241    \tab 12.102    \tab 11.310 \cr
6 \tab 100001\_at     \tab 5.803     \tab 6.148     \tab 5.969 \cr
}

}
                                             
\references{
  Irizarry RA, Hobbs B, Collin F, Beazer-Barclay YD,
  Antonellis KJ, Scherf U, Speed TP. Exploration, normalization, and
  summaries of high density oligonucleotide array probe level data. 
  Biostatistics. 2003 Apr;4(2):249-64.}

\note{
  The processing time of this function, as well as the size of the output value, 
  increases linearly with the number of .CEL files supplied.
}

\author{ 
Rudi Alberts <r.alberts@rug.nl>, 
Ritsert Jansen <r.c.jansen@rug.nl>
}

\seealso{
  \code{\link{rma.preprocessing.plussummary}}
}

\examples{
\dontrun{

# We set the directory where the .CEL files are placed
directory <- "C:/data/my_project/celfiles"

# We change the working directory to the path we want our output file 
#  to be placed
setwd("C:/data/my_project/celfiles")

# We pre-process our .CEL files 
celfiles <- c("bxd1a.cel","bxd2a.cel","bxd5a.cel")    
probesignals <- rma.preprocessing( directory, celfiles,
    "probesignals.csv" )

# We check the size and headers of the output
dim(probesignals)               
head(probesignals)

}
}
