\name{rma.preprocessing.plussummary}
\alias{rma.preprocessing.plussummary}
\title{ Preprocessing Affymetrix .CEL files using RMA including summary }

\description{
  Applies the RMA preprocessing method to Affymetrix .CEL files.
  Output is obtained on probeset level (with summary).
}

\usage{
rma.preprocessing.plussummary( directory, celfiles, filename,
                               cdffile, cdfpackage )
}

\arguments{

  \item{directory}{ path of the directory where files are located.\cr
  In this function, input directory where .CEL files are located. }
  
  \item{celfiles}{ name of the .CEL files that should be processed. }
  
  \item{filename}{ (optional) path of the file where the data will be loaded from or save to.\cr
  In this function, output file where results are stored (*.txt = tab
    separated values). Default is NULL. }
    
  \item{cdffile}{ (optional) filename argument with a custom Affymetrix chip 
                   description file (CDF). Default is NULL.}
  
  \item{cdfpackage}{ (optional) filename argument with a package containing the CDF file name. 
                   Default is NULL.}
}

\details{
  \code{\link{rma.preprocessing.plussummary}} does RMA (Robust Multichip Analysis)
  background correction, quantile normalization and RMA summarization. 
  Probe set level results are stored.
}

\value{
  Returns a \eqn{Biobase Expression Set (exprSet)} (and creates a .txt file) 
  where on each line the summarized probe set values for all .CEL files supplied are given. \cr
  Example with 3 .CEL files:
  \preformatted{
> ps.summ <- rma.preprocessing.plussummary( directory, "ps.summ.txt" )
> class(probeset.summary)
[1] "exprSet"
attr(,"package")
[1] "Biobase"
> probeset.summary
Expression Set (exprSet) with 
        12488 genes
        3 samples
        phenoData object with 1 variables and 3 cases
        varLabels
                sample: arbitrary numbering
 
}
  Beginning of the .txt output file: \cr

  \tabular{lrrr}{
  \tab bxd1a.cel	\tab bxd2a.cel	\tab bxd5a.cel \cr
100001\_at	\tab 8.9134282866486	\tab 9.04991691853465	\tab 8.89581017824296 \cr
100002\_at	\tab 8.57571818513605	\tab 8.67811529853219	\tab 8.6293385031351 \cr
100003\_at	\tab 7.58387242256224	\tab 7.68779109466728	\tab 7.59685409369385 \cr
100004\_at	\tab 9.74602889100652	\tab 9.74602889100652	\tab 9.74602889100652 \cr
100005\_at	\tab 10.0287767421828	\tab 10.0367757465467	\tab 9.98787086359905 \cr
}
 
}

\references{
  Irizarry RA, Hobbs B, Collin F, Beazer-Barclay YD,
  Antonellis KJ, Scherf U, Speed TP. Exploration, normalization, and
  summaries of high density oligonucleotide array probe level data. Biostatistics. 2003 Apr;4(2):249-64.}

\note{
  The processing time of this function, as well as the size of the output value, 
  increases linearly with the number of .CEL files supplied.
}

\author{ 
Rudi Alberts <r.alberts@rug.nl>, 
Ritsert Jansen <r.c.jansen@rug.nl>
}
\seealso{
  \code{\link{rma.preprocessing}}
}
\examples{
\dontrun{

# We set the directory where the .CEL files are placed
directory <- "C:/data/my_project/celfiles"

# We change the working directory to the path we want our output file to be placed
setwd("C:/data/my_project/celfiles")

# We pre-process our .CEL files 
ps.summ <- rma.preprocessing.plussummary( directory, "ps.summ.txt" )

# We check the content of the output
ps.summ     

}
}
