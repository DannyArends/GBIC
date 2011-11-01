\name{experimentDesignTable}
\alias{experimentDesignTable}
\title{ Make experiment table based design matrix }
\description{
  This function generates two \code{.csv} files which descibe how samples are allocated
   samples into different conditions and paired on arrays. 
}
\usage{
  experimentDesignTable( array.allocation, condition.allocation, 
                         nEnvFactors, nLevels, Level, fileName,envFactorNames,
                         directory )
}
\arguments{
   \item{array.allocation}{ 
      matrix with nArray rows and nRIL columns.
      Elements of 1/0 indicate this RIL is/not selected 
      for this array. 
  }
  \item{condition.allocation}{ 
      matrix with nCondition rows and nRIL columns.
      Elements of 1/0 indicate this RIL is/not selected 
      for this condition. 
  }
   \item{nEnvFactors}{
      number of environmental factors, an integer.
      When \code{nEnvFactors} is between 0 and 3, all main factor(s) and
      interacting facotr(s) will be included.
  }
  \item{nLevels}{
      number of levels for each factor, a vector with each
      component being integer. The length of it should equal
      \code{nEnvFactors}.
  }
 \item{Level}{
      a list which specifies the levels for each factor in the
      experiment. There are in total \code{nEnvFactors} elements in the list and each
      element correpsond to certain envrironmental factor. The emlemet is a
      vector describing all levels of the environmental factor.
      default setting for the level of each factor is 1, 2, ...nLevels[i]. (Here
      nLevels[i] is the \emph{i}th element of nLevels, which gives the total number of
      levels for \emph{i} environmental facotor).
  }
  \item{fileName}{
         the final optimal design table(s) in \code{csv} format and a plot (in \code{png} 
         format) of the all scores during SA process (if \code{plotScores} = T) will be 
         produced. The users can specify the table and plot name by setting 
         \code{fileName}.
         If \code{NULL} (default) it produces files starting with \code{"myDesignGG"}.
  }
  \item{envFactorNames}{
        a vector with names for all environmental factor(s). For example, for the
        experiment with two environmental factors of temperature and drug treatment:
        \code{envFactorNames <- c( "Temperature", "Dosage" ) } \cr
        Default = \code{NULL}, then the output will use "F1" and "F2" to indicate the
        environmental factors.
  }
  \item{directory}{
       It tells where the resulting optimal design tables are to be stored.
       If \code{NULL} (default), it will use the currect working directory.
  } 
}
\details{
 Based on \code{nEnvFactors} and \code{nLevels}, \code{nConditions} is calculated. 
}
\value{
  Two tables report the results: table "pair design" which is only used for 
  two-channel experiments and describes how samples are paired together on the 
  slide (e.g. microarray chip), and table "envi-ronment design" which is used 
  when there are more environments evolved in the experiment. With these two 
  tables, the experimenters can set up the environmental treatment and 
  follow-up profiling measurement.  \cr
  Examples: \cr
  1. conditionDesign.csv \cr
\tabular{lccccccc}{
            \tab Temperature  \tab Cell Type   \tab \tab \tab Selected Samples \tab \tab   \cr 
condition1  \tab 15 \tab A \tab RIL28 \tab RIL81 \tab RIL18 \tab RIL61 \tab \cr  
condition2  \tab 24 \tab A \tab RIL72 \tab RIL40 \tab RIL83 \tab RIL44 \tab RIL10 \cr 
condition3  \tab 29 \tab A \tab RIL22 \tab RIL89 \tab RIL3 \tab RIL30 \tab RIL58 \cr 
condition4  \tab 15 \tab B \tab RIL70 \tab RIL47 \tab RIL4 \tab RIL59 \tab \cr  
condition5  \tab 24 \tab B \tab RIL93 \tab RIL97 \tab RIL49 \tab RIL14 \tab \cr
} 
  2. arrayDesign.csv  \cr
 \tabular{lrr}{
            \tab Channel 1 \tab Channel 2 \cr 
array1      \tab RIL28    \tab RIL92 \cr    
array2      \tab RIL70    \tab RIL47 \cr    
array3      \tab RIL22    \tab RIL89 \cr    
array4      \tab RIL45    \tab RIL15 \cr    
array5      \tab RIL52    \tab RIL41 \cr
}
}
\references{ 
  Y. Li, R. Breitling and R.C. Jansen. Generalizing genetical
       genomics: the added value from environmental perturbation, Trends Genet
       (2008) 24:518-524. \cr
  Y. Li, M. Swertz, G. Vera, J. Fu, R. Breitling, and R.C. Jansen. designGG:
  An R-package and Web tool for the optimal design of genetical genomics
  experiments. (submitted) \cr
  http://gbic.biol.rug.nl/designGG

}
\author{  Yang Li <yang.li@rug.nl>, Gonzalo Vera <gonzalo.vera.rodriguez@gmail.com> \cr
  Rainer Breitling <r.breitling@rug.nl>, Ritsert Jansen <r.c.jansen@rug.nl> 
 }
\note{ 
    The optimal design results are described in two tables. One is called "array
    design" which is only used for two-channel experiments. It describes how 
    samples are paired together on the slide (e.g. microarray chip). The other 
    table is  called "condition design" which is used when there is more than
    one environmental factor involved in the experiment. Each cell in condition 
    design table represents a combination of different levels of environmental
    factors and the selected sample names (e.g. RIL names) for this condition
    are shown. Based on these two tables, the experimenters can set up the 
    environmental treatment and follow-up profiling measurement.
}
\seealso{  \code{\link{designGG}}, \code{\link{exampleArrayDesignTable}}, \code{\link{exampleConditionDesignTable}} }
\keyword{method}
