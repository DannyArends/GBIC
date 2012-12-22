\name{splitProbeSetTraits}
\alias{splitProbeSetTraits}
\title{ Split Affymetrix signals into files per probe set}
\description{
  After preprocessing .CEL files with \code{\link{rma.preprocessing}}
  the result is a \code{data.frame} with the signals (traits) obtained from 
  several individuals.
  This function helps to create, per probe set, one file with the trait values of each probe. 
  Later, these files can be loaded and used to analyze only the probe set of interest.
}
\usage{
splitProbeSetTraits(affy.signals, directory)
}
\arguments{
  \item{affy.signals}{ a \code{data.frame} where on each line the probe signal 
                      values for all .CEL files supplied are given. \cr
                      See also return value of \code{\link{rma.preprocessing}}.
  } 
  \item{directory}{ path of the directory where new csv files will be located. }
}
\details{
  The filename of the files created at the \code{directory} indicated 
  have the following format:\cr
  \code{directory/traits.probeset_name.csv}\cr
  
  As it is possible that the probe set name contains illegal characters, the
  \eqn{candidate} file names are filtered and modified if needed with the function
  \code{\link{auxCleanFilename}}.
}

\value{
This function has no return value. \cr
It creates as many files as probe sets are provided 
in \code{affySignals}. The following example illustrates how the files are
created and which changes are introduced:\cr
The next rows on \code{affySignals}:
\preformatted{
affySignals[ affySignals\$probeset == "91241/01_at", ]
}
\tabular{llrrrr}{
      \tab probeset     \tab bxd1a.cel  \tab bxd1b.cel \tab bxd2a.cel \tab bxd2b.cel \cr
98412 \tab 91241/01\_at \tab     7.834  \tab   7.968   \tab     8.056 \tab     8.517 \cr
98413 \tab 91241/01\_at \tab     8.237  \tab   8.233   \tab     8.796 \tab     8.869 \cr
98414 \tab 91241/01\_at \tab    11.353  \tab   11.273  \tab    11.360 \tab    11.526 \cr
98415 \tab 91241/01\_at \tab    11.644  \tab   11.044  \tab    10.969 \tab    11.137 \cr
98416 \tab 91241/01\_at \tab     7.437  \tab    7.049  \tab     7.486 \tab     6.012 \cr
}
generates a file called: traits.\eqn{91241\_01\_at}{91241_01_at}.csv, with the following contents:
\tabular{lrrrr}{
              \tab  BXD1A  \tab  BXD1B \tab  BXD2A \tab  BXD2B \cr
91241\_01\_at1 \tab  7.834  \tab  7.968 \tab  8.056 \tab  8.517 \cr
91241\_01\_at2 \tab  8.237  \tab  8.233 \tab  8.796 \tab  8.869 \cr
91241\_01\_at3 \tab 11.353  \tab 11.273 \tab 11.360 \tab 11.526 \cr
91241\_01\_at4 \tab 11.644  \tab 11.044 \tab 10.969 \tab 11.137 \cr
91241\_01\_at5 \tab  7.437  \tab  7.049 \tab  7.486 \tab  6.012 \cr
}
}

\author{ 
Gonzalo Vera <g.vera.rodriguez@rug.nl>
}
\note{ 
  Usually, due to the high amount of probe sets obtained from Affymetrix .CEL files,
  a high amount of files are generated in a single directory. Depending on the 
  OS and the machine used, this can result in longer waiting times while the OS 
  retrieves the complete list of files. Therefore, we advise against processing at 
  one time all the probe sets.\cr
  Splitting the input \code{affySignals} in a few thousands probe sets is recommended,
  mainly when working with Windows environments. \cr
  \code{directory} must end without "/".
}
\seealso{
  \code{\link{rma.preprocessing}}, \code{\link{auxCleanFilename}}
}
\examples{
\dontrun{

# We get the signals from .CEL files
directory <- "C:/data/my_project/celfiles"
affySignals <- rma.preprocessing( directory )

# We split it in small files x trait
directory <- "C:/data/my_project/traitfiles"
splitProbeSetTraits( affySignals, directory )

}
}
