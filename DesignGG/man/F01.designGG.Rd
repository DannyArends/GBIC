\name{designGG}
\alias{designGG}
\alias{design genetical genomics experiment}
\title{Optimal design for genetical genomics experiments }
\description{
  Main function to search and display A- and D- optimal designs for single or
  dual channel genetical genomics experiments. Simulated annealing or Metropolis 
  Hastings is used to find the best design.
}
\usage{
  designGG <- function( genotype, nSlides, nTuple, 
			nEnvFactors, nLevels, Level=NULL, 
			bTwoColorArray=T, initial=NULL, 
			weight=1, region=NULL, optimality="A", 
			method="SA", nIterations=3000, n.search=2, 
			endTemp=1e-10, startTemp=1, maxTempStep=0.9,
                     	plotScores=T, directory=NULL, fileName=NULL, 
                     	envFactorNames=NULL, writingProcess=F )
                     
}
\arguments{
  \item{genotype}{
      genotype data: a matrix nMarker x nRILs with two allels being 0 and 1.
  }
  \item{nSlides}{
      total number of slides available for experiment
  }
  \item{nTuple}{
      average number of RILs put into each condition nTuple > 1, a real number.
      if nTuple < 1, stop,  \cr
      \code{warning message:"The number slides is too less to perform the experiment."}
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
      experiment. There are in total nEnvFactors elements in the list and each 
      element correpsond to certain envrironmental factor. The emlemet is a 
      vector describing all levels of the environmental factor.     
      default setting for the level of each factor is 1,2,....nlevels[i]. (Here
      nLevels[i] is the ith element of nLevels, which tells the total number of 
      levels for ith environmental facotor) 
  }
  \item{bTwoColorArray}{
      binary variable indicating experiment type: \cr
      bTwoColorArray = T, dual channel experiment \cr
      bTwoColorArray = F, single channel experiment 
  }
  \item{initial}{
        The starting design matrix for the algorithm. If specified this should
        be a list with 2 matrices: \cr
        condition.allocation: allocate RILs into different conditional
        (nrow = nCondition, ncol= nRILs) \cr
        array.allocation: pair RILs into sldies
        (nrow = nSlide, ncol = nRILs)       \cr
        However, the algorithm does not require that a starting matrix is
        specified. Default = NULL.
  }
  \item{weight}{
      a vector with length of variablNumber which is calculated
      from function variablNumber, default = 1 which means the parameters to be 
      estimated are equally important during optimization.
  }
  \item{region}{
      genome region of biological interest. \cr
      default = NULL (means all genome will considered).
  }
  \item{optimality}{
      type of optimality, i.e. "A" (A-optimality) or "D"
      (D-optimality). A-optimality minimizes $Trace((X'X)^{-1})$,
      which corresponds to minimum average variance of the
      parameter estimates. D-optimality minimizes $det(X'X)^{-1}$,
      which corresponds to minimum generalized variance of the
      parameter estimates.
  }
  \item{method}{
      Method for searching for an optimal design. "SA" uses simulated annealing.
      "MH" uses Metropolis Hasting,   default: "SA".
  }
  \item{nIterations}{
        Number of iterations of the simulated annealing method. default = 3000.
  }
  \item{n.search}{
        number of times for simulated annealing optimaization with different
        initial design , default = 2. 
        Here it is suggested to be between 1 and 5. 
        It should not to be too large because of computational burden.
  }
  \item{endTemp}{
        ending temperature of simulated annealing process.
        An important optimization parameter.
        Setting this parameter closer to zero. defaut = $1e^{-10}$.
  }
  \item{startTemp}{
        starting temperature of simulated annealing process.  default = 1.
  }
  \item{maxTempStep}{
        maximum  temperature decreasing step for simulated annealing process.
        The parameter ensures that the multiplicative cooling factor is not
        smaller than that. If nIterations is too small, the preferred final
        temperature (endTemp) may not be reached.See Wit and McClure (2004)
        for details.Here default = 0.9.
  }
  \item{plotScores}{
         If TRUE (the default) it produces a plot of the design found by SA
         using the function plotAllScores.
  }
  \item{directory}{
       It tells where the resulting optimal design tables are to be stored.
       If NULL (default), it will take currect working directory.
  } 
  \item{fileName}{
         The final optimal design table(s) in csv format and a plot(in png 
         format) of the all scores during SA process (if plotScores=T) will be 
         produced. The users can specify the table and plot name by setting 
         fileName.
         If NULL (the default) it produces files starting with "myDesignGG".
  }
  \item{envFactorNames}{
        a vector with names for all environmental factor(s). For example, for the
        experiment with two environmental factors of temperature and cell types
        envFactorNames <- c( "Temperature", "Dosage" )
        default = NULL, then the output will use "var1" and "var2" to indicate the
        environmental factors.
  }
  \item{writingProcess}{
         If TRUE, it prints how much computation work has been finished in a
         file called "processing.txt".
         defaut=FALSE
  }
}

\details{
  Given the genetic information of samples available for the experiment
  (genotype) and the information about experimental settings (nEnvFactors,
  nSlides, nLevels etc.), the algorithm searches foran A-optimal or D-optimal 
  (see "optimality") using simulated annealing (see "method"). A plot of
  the scores at each iterations can also be given using the plotAllScores 
  function. It also contains a number of the arguments: "region" to specify the 
  genome region that are major interest to experimenters, "weight" to define
  the weight of  genetic, environmental factor and interaction terms.
  Example parameter settings:
  Suppose to design an experiment with two environmental factors (F1, F2) and 
  there are two diffferent levels for each environment. The levels are 16
  and 24 for F1, and 5 anf 10 for F2. Thus the following command can be used:\cr
      > nEnvFactors 	<- 2       \cr
      > nLevels      	<- c ( 2, 2 ) \cr
      > levels 	      <- list ( c(16, 24), c(5, 10) )  \cr
  The length of parameter "weight" is dependent on the number of environmental
  factors: \cr
  When nEnvFactor = 0, \cr
    weight is 1 as there is only one parameter of interest (genotype).\cr
  When nEnvFactor = 1, \cr 
    weight = c( $w_{Q}$, $w_{F1}$, $w_{QF1}$ ) \cr
  When nEnvFactor = 2, \cr 
    weight = c( $w_{Q}$, $w_{F1}$, $w_{F2}$, $w_{QF1}$, $w_{QF2}$, $w_{F1F2}$, $w_{QF1F2}$) \cr
  When nEnvFactor = 3, \cr
    weight = c( $w_{Q}$, $w_{F1}$, $w_{F2}$, $w_{F2}$,
                $w_{QF1}$, $w_{QF2}$, $w_{QF3}$, $w_{F1F2}$, $w_{F1F3}$, $w_{F2F3}$,
                $w_{QF1F2}$, $w_{QF1F3}$, $w_{QF2F3}$, $w_{QF1F2F3}$ ) \cr    
  Here $w_{Q}$ represent the weight for genotype effect,  $w_{F1}$ represent the 
  weight for F1 effect and $w_{QF1}$ represent the weight for interaction between
  genotype and F1 effect, ect. 
}
\value{

  An array design table (arrayDesign.csv) and a condition design table (
    conditionDesign.csv) are generated.

}
\references{
     Y. Li, R. Breitling and R.C. Jansen "Generalizing genetical 
     genomics: the added value from environmental perturbation", Trends Genet
     (2008) 24:518-524.
     
     E. Wit and J. McClure  "Statistics for Microarrays: Design, Analysis
      and Inference" (2004) Chichester: Wiley.
}
\author{  
 Yang Li <yang.li@rug.nl>, Gonzalo Vera <g.vera.rodriguez@rug.nl> \cr
 Rainer Breitling <r.breitling@rug.nl>, Ritsert Jansen <r.c.jansen@rug.nl> 
}

\seealso{ 
  \code{\link{initialDesign}},
  \code{\link{designScore}}, 
  \code{\link{updateDesign}}, 
  \code{\link{acceptanceProbability}}, \cr
  \code{\link{experimentDesignTable}}, 
  \code{\link{plotAllScores}}
}
\examples{
  library(designGG)
  #load genotype
  data(genotype)
  #Example:  single-channel experiment with 2 envrionmental factors,
  #each with 5 level, and there will be four samples per condition(nTuple=4).
  optimalDesign <- designGG ( genotype, nSlides=NULL, nTuple=4, nEnvFactors=2,
                        nLevels=c(5,5),Level=NULL,  bTwoColorArray=F,
                        initial=NULL, weight=1, region=NULL, optimality="A", 
                        method="SA", nIterations=100, n.search=2, endTemp=1e-10,
                        startTemp=1, maxTempStep=0.9, plotScores=T, 
                        directory=NULL, fileName=NULL, envFactorNames=NULL, 
                        writingProcess=F )
  #example 2:  dual-channel experiment with 2 envrionmental factors,
  #each with 5 level. There are 50 slides available.
  optimalDesign <- designGG ( genotype, nSlides=50, nTuple=NULL, nEnvFactors=2,
                        nLevels=c(5,5),Level=NULL,  bTwoColorArray=T,
                        initial=NULL, weight=1, region=NULL, optimality="A", 
                        method="SA", nIterations=100, n.search=2, endTemp=1e-10,
                        startTemp=1, maxTempStep=0.9, plotScores=T, 
                        directory=NULL, fileName=NULL, envFactorNames=NULL, 
                        writingProcess=F )
  #result
  optimalDesign$arrayDesign
  optimalDesign$conditionDesign
  plotAllScores(optimalDesign$plot.obj)
}
\keyword{ method }
