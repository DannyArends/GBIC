# Function name: rma.preprocessing
# Author: Rudi Alberts <r.alberts@rug.nl>
# Maintainer: Gonzalo Vera <g.vera.rodriguez@rug.nl>
# Version: 1.1.1
# Date: 02 July 2007

rma.preprocessing <- function(directory, celfiles, filename=NULL, cdffile = NULL, cdfpackage = NULL){
  library(affy)

  # readin the specified 'celfiles' in directory
  A   <-  ReadAffy( celfile.path=directory, filenames=celfiles ) 
 
  # If supplied, we load the custom CDF files
  if (!is.null(cdffile) & !is.null(cdfpackage)){
      library(cdfpackage)
      A@cdfName<-'cdffile'
  }
  
  # RMA background correction
  A2  <-  bg.correct.rma( A )
  # RMA quantile normalization
  A3  <-  normalize.AffyBatch.quantiles( A2 )  
  
  # take log2 and add probe names
  pmperprobe  <- cbind(probeset=probeNames(A3), round( log2(data.frame(pm(A3))),3 ) )

  if (!is.null(filename)){ # write output to file
    write.table(pmperprobe, filename, sep=",", row.names=T)
  }
  pmperprobe
}

