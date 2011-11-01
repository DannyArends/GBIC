# Originally made to assist splitTraits func
# Function name: auxCleanFilename
# Author: Gonzalo Vera <g.vera.rodriguez@rug.nl>
# Maintainer: Gonzalo Vera <g.vera.rodriguez@rug.nl>
# Version: 1.0.0
# Date: 1 May. 2007
# TODO: Implement something more efficient with string2list and apply

# input: a string with a candidate filename
# output: a cleaned filename where illegal carachters are replaced by valid ones
# Example:
# [1] "? [ ] / \\ = + < > : ; \" ,"
# [1] "# { } _ _ - _ { } _ _ ' _"

# Check that the name doesn't contains reserved caracthers like ? [ ] / \ = + < > : ; " ,
auxCleanFilename <- function ( filename=NULL )
{
  if( is.null(filename) )
  {
      stop(" In gg.auxCleanFilename: You must supply a filename" )
  }
  invalidChar <- c( "\\?", "\\[", "]", "/", "\\\\", "=", "\\+", "<", ">", ":", ";", "\"", "," )
  replaceChar <- c( "#",   "{",   "}", "_", "_", "-", "_", "{", "}", "_", "_", "'",  "_" )

  cleanFilename <- filename
  for( idx in 1:length(invalidChar) )
  {
    cleanFilename <- gsub( invalidChar[ idx ],
                           replaceChar[ idx ], cleanFilename )
  }      
  
  cleanFilename
}