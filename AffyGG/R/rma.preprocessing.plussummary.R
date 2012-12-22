# Function name: rma.preprocessing.plussummary
# Author: Rudi Alberts <r.alberts@rug.nl>
# Maintainer: Gonzalo Vera <g.vera.rodriguez@rug.nl>
# Version: 1.1.0
# Date: 1 Feb. 2006


rma.preprocessing.plussummary <- function( directory, celfiles, filename=NULL, 
                                           cdffile = NULL, cdfpackage = NULL )
{
  library( affy )

  # readin the specified 'celfiles' in directory
  A   <-  ReadAffy( celfile.path=directory, filenames=celfiles )
  
  # If supplied, we load the custom CDF files
  if (!is.null(cdffile) & !is.null(cdfpackage)) 
  {
      library(cdfpackage)
      A@cdfName<-'cdffile'                         
  }
  
  # RMA background correction, normalization and summary
  eset  <- rma(A)             
  
  if (!is.null(filename)) 
  {                                   # write output to file
    write.exprs( eset, file=filename) # eset belongs to class exprSet
  }
  eset
}
