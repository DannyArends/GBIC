\name{experimentDesignTable}
\alias{experimentDesignTable}
\title{ Make experiment table based design matrix }
\description{
  This function generates two .CSV files which descibe how samples are allocated
   samples into different conditions and paired on arrays. 
}
\usage{
  experimentDesignTable( array.allocation, condition.allocation, 
                         nEnvFactors, nLevels, Level, fileName,envFactorNames,
                         directory )
}
\arguments{
  \item{array.allocation}{ 
      matrix iwth nArray rows and nRIL columns.
      elements of 1/0 indicates this RIL is/not selected 
      for this array. 
  }
  \item{condition.allocation}{ 
      matrix iwth nCondition rows and nRIL columns.
      elements of 1/0 indicates this RIL is/not selected 
      for this condition. 
  }
  \item{nEnvFactors}{ 
      number of environmental factors, integer.
      When nEnvFactors is bewteen 0 and 3, all main factor and
      interacting facotrswill be included.
  }
  \item{nLevels}{ 
      number of levels for each factor, a vector with each 
      component being integer. The length of it should equal 
      nEnvFactors.
  }
  \item{Level}{ 
      a list which specifies the levels for each factor in the 
      experiment. The length of it should be nEnvFactors.
      default: Level[i]=(1,2,3..nLevels[i])
  }
  \item{fileName}{
         The final optimal design table(s) in csv format  will be produced. 
         The users can specify the table name by setting fileName.
         If NULL (the default) it produces files starting with "myDesignGG".
  }
  \item{envFactorNames}{
        a vector with names for all environmental factor(s). For example, for the
        experiment with two environmental factors of temperature and cell types
        envFactorNames <- c( "Temperature", "CellType" )
        default=NULL, then the output will use "var1" and "var2" to indicate the
        environmental factors.
  }
  \item{directory}{
       It tells the function where the resulting optimal design tables are to 
       be stored. 
       If NULL (default), the currect working directory is taken.
  } 
}
\details{
 Based on nEnvFactors and nLevels, nConditions is calculated. 
}
\value{
  Two tables report the results: table "pair design" which is only used for 
  two-channel experiments and describes how samples are paired together on the 
  slide (e.g. microarray chip), and table "envi-ronment design" which is used 
  when there are more environments evolved in the experiment. With these two 
  tables, the experimenters can set up the environmental treatment and 
  follow-up profiling measurement.  \cr
  Exmaples: \cr
  1. conditionDesign.csv \cr
 \tabular{lrrrrrrr}{
            \tab Temperature  \tab Cell Type   \tab \tab \tab \tab \tab   \cr 
condition1  \tab 15 \tab A \tab RIL28 \tab RIL81 \tab RIL18 \tab RIL61 \tab \cr  
condition2  \tab 24 \tab A \tab RIL72 \tab RIL40 \tab RIL83 \tab RIL44 \tab RIL10 \cr 
condition3  \tab 29 \tab A \tab RIL22 \tab RIL89 \tab RIL3 \tab RIL30 \tab RIL58 \cr 
condition4  \tab 15 \tab B \tab RIL70 \tab RIL47 \tab RIL4 \tab RIL59 \tab \cr  
condition5  \tab 24 \tab B \tab RIL93 \tab RIL97 \tab RIL49 \tab RIL14 \tab \cr
} 
  2. arrayDesign.csv  \cr
 \tabular{lrr}{

            \tab Channel1 \tab Channel2 \cr 
array1      \tab RIL28    \tab RIL92 \cr    
array2      \tab RIL70    \tab RIL47 \cr    
array3      \tab RIL22    \tab RIL89 \cr    
array4      \tab RIL45    \tab RIL15 \cr    
array5      \tab RIL52    \tab RIL41 \cr
}
}
\references{ 
   Y. Li, R. Breitling and R.C. Jansen "Generalizing genetical 
     genomics: the added value from environmental perturbation", Trends Genet
     (2008) 24:518-524. 
}
\author{  
 Yang Li <yang.li@rug.nl>, Gonzalo Vera <g.vera.rodriguez@rug.nl> \cr
 Rainer Breitling <r.breitling@rug.nl>, Ritsert Jansen <r.c.jansen@rug.nl> 
}
\note{ 
    The optimal design results are described in two tables. One is called "array
    design" which is only used for two-channel experiments. It describes how 
    samples are paired together on the slide (e.g. microarray chip). The other 
    table is  called "condition design" which is used when there are more than
    one environmental factors evolved in the experiment. Each cell in condition 
    design table represents a combination of different levels of environmental
    factors and the selected sample names (e.g. RIL names) for this condition
    are shown. Based on these two tables, the experimenters can set up the 
    environmental treatment and follow-up profiling measurement.
}
\seealso{ See Also \code{\link{designGG}},\code{\link{exampleArrayDesignTable}},
\code{\link{exampleConditionDesignTable}} }
\keyword{method}
