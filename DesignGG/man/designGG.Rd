\name{designGG}
\alias{designGG}
\alias{design genetical genomics experiment}
\title{Optimal design for genetical genomics experiments }
\description{
  Main function to search and display A- and D- optimal designs for single- or
  two-channel genetical genomics experiments. Simulated annealing or Metropolis 
  Hastings used to find the best design.
}
\usage{
  designGG(genotype, nSlides, nTuple, nEnvFactors, nLevels,
           Level=NULL, bTwoColorArray=TRUE, initial=NULL, weight=1,
           region=NULL, optimality="A", method="SA", nIterations=3000,
           n.search=2, endTemp=1e-10, startTemp=1, maxTempStep=0.9,
           plotScores=TRUE, directory=NULL, fileName=NULL, 
           envFactorNames=NULL, writingProcess=FALSE )                     
}
\arguments{
  \item{genotype}{
      genotype data: a nMarker-by-nRILs matrix with two allels being 0 and 1 (or A and B)
      or three allels being 0, 0.5 and 1 (or, A, H, and B), where 0.5 (or H) represents
      heterozygous allele.
  }
  \item{nSlides}{
      total number of slides available for the experiment.
  }
   \item{nTuple}{
      average number of RILs to be assigned onto each condition. \cr
      \code{nTuple} should be a real number which is larger than 1.\cr
      If \code{nTuple} < 1, the algorithm will stop and show the message, \cr
      \code{warning: "The number of slides is too small to perform the experiment."}
  }
  \item{nEnvFactors}{
      number of environmental factors, integer.
      When \code{nEnvFactors} is between 0 and 3, all main factor(s) and
      interacting facotr(s) will be included.
  }
  \item{nLevels}{
      number of levels for each factor, a vector with each
      component being an integer. The length of it should equal
      \code{nEnvFactors}.
  }
 \item{Level}{
      a list which specifies the levels for each factor in the
      experiment. There are in total \code{nEnvFactors} elements in the list and each
      element corresponds to a certain environmental factor. The element is a
      vector describing all levels of the environmental factor.
      default setting for the level of each factor is 1, 2, ...., nlevels[i]. (Here
      nLevels[i] is the \emph{i}th element of nLevels, which gives the total number of
      levels for \emph{i} environmental factor).
  }
  \item{bTwoColorArray}{
      binary variable indicating experiment type: \cr
      \code{bTwoColorArray <- TRUE} \#for dual channel experiment \cr
      \code{bTwoColorArray <- FALSE} \#for single channel experiment 
  }
  \item{initial}{
        the starting design matrix for the algorithm. If specified, this should
        be a list with 2 matrices: \cr
        \code{condition.allocation}: allocate RILs into different conditional
        (nrow = nCondition, ncol= nRILs) \cr
        \code{array.allocation}: pair RILs into sldies
        (nrow = nSlide, ncol = nRILs)       \cr
        However, the algorithm does not require that a starting matrix is
        specified. Default = \code{NULL}.
  }
  \item{weight}{
      a vector with length of \code{variableNumber} which is calculated
      from function \code{variableNumber}. Default = 1 (which means the parameters to be 
      estimated are all equally important during optimization). See details below.
  }
  \item{region}{
      genome region of biological interest. Default = \code{NULL} (which means
      the entire genome will considered).
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
      method for searching for an optimal design. "SA" uses simulated annealing.
      "MH" uses Metropolis Hasting. Default = "SA".
  }
  \item{nIterations}{
        number of iterations of the simulated annealing method. Default = 3000.
  }
  \item{n.search}{
        number of times for simulated annealing optimaization with different
        initial design, default = 2. 
        Here it is suggested to be between 1 and 5. 
        It should not to be too large because of the reaching computational burden.
  }
  \item{endTemp}{
        ending temperature of simulated annealing process.
        An important optimization parameter.
        Default = $1e^{-10}$.
  }
  \item{startTemp}{
        starting temperature of simulated annealing process.  Default = 1.
  }
  \item{maxTempStep}{
        maximum temperature decreasing step for simulated annealing process.
        The parameter ensures that the multiplicative cooling factor is not
        smaller than that. If \code{nIterations} is too small, the preferred final
        temperature (endTemp) may not be reached. See Wit and McClure (2004)
        for details. Default = 0.9.
  }
  \item{plotScores}{
         If \code{TRUE} (default) it produces a plot of the optimazation by SA
         using the function plotAllScores.
  }
  \item{directory}{
       It tells where the resulting optimal design tables are to be stored.
       If \code{NULL} (default), it will take currect working directory.
  } 
  \item{fileName}{
         the final optimal design table(s) in \code{csv} format and a plot (in \code{png} 
         format) of all scores during SA process (if \code{plotScores} = T) will be 
         produced. The users can specify the table and plot name by setting 
         \code{fileName}.
         If \code{NULL} (default) it produces files starting with \code{"myDesignGG"}.
  }
 \item{envFactorNames}{
        a vector with names for all environmental factor(s). For example, for the
        experiment with two environmental factors of temperature and drug treatment:
        \code{envFactorNames <- c( "Temperature", "Dosage" ) } \cr
        Default = \code{NULL}, then the output will use \code{"F1"} and \code{"F2"} to indicate the
        environmental factors.
  }
  \item{writingProcess}{
         If TRUE, it prints how much computation work has been finished in a
         file called \code{"processing.txt"}.
         Default = \code{FALSE}.
  }
}

\details{
  Given the genetic information of samples available for the experiment
  (genotype) and the information about experimental settings (\code{nEnvFactors},
  \code{nSlides},\code{nLevels} etc.), the algorithm searches for an A-optimal or D-optimal 
  (see \code{optimality}) using simulated annealing (see \code{method}). A plot of
  the scores at each iterations can also be given using the \code{plotAllScores} 
  function. \cr
  It also contains a number of the arguments:\cr
  \code{region} is used to specify the 
  genome region that are of major interest to experimenters. \cr
  \code{weight} is used to define
  the weight of genetic and environmental factors,  and interaction terms. Prior 
  knowledge about expected effect sizes of interesting factors can also be 
  incorporated as \code{weight} parameters for the algorithm. The weight is 
  inversely proportional to the expected effect size of the corresponding parameter.
  Example parameter settings:
  Suppose to design an experiment with two environmental factors (F1, F2) and 
  there are two diffferent levels for each environment. The levels are 16
  and 24 for F1, and 5 and 10 for F2. Thus the following command can be used:\cr
  \code{ nEnvFactors 	<- 2 }      \cr
  \code{ nLevels     	<- c ( 2, 2 )} \cr
  \code{ levels 	    <- list ( c(16, 24), c(5, 10) ) } \cr
  The length of parameter \code{weight} is dependent on the number of environmental
  factors: \cr
  When \code{nEnvFactor} = 0, \cr
    \code{weight} is 1 as there is only one parameter of interest (genotype).\cr
  When \code{nEnvFactor} = 1, \cr 
    \code{weight} = c( $w_{Q}$, $w_{F1}$, $w_{QF1}$ ) \cr
  When \code{nEnvFactor} = 2, \cr 
    \code{weight} = c( $w_{Q}$, $w_{F1}$, $w_{F2}$, $w_{QF1}$, $w_{QF2}$, $w_{F1F2}$, $w_{QF1F2}$) \cr
  When \code{nEnvFactor} = 3, \cr
    \code{weight} = c( $w_{Q}$, $w_{F1}$, $w_{F2}$, $w_{F2}$,
                $w_{QF1}$, $w_{QF2}$, $w_{QF3}$, $w_{F1F2}$, $w_{F1F3}$, $w_{F2F3}$,
                $w_{QF1F2}$, $w_{QF1F3}$, $w_{QF2F3}$, $w_{QF1F2F3}$ ) \cr    
  Here $w_{Q}$ represents the weight for genotype effect,  $w_{F1}$ represent the 
  weight for F1 effect and $w_{QF1}$ represent the weight for interaction between
  genotype and F1 effect, etc. 
}
\value{

  An array design table (arrayDesign.csv) and a condition design table (
    conditionDesign.csv) are generated.

}
\references{
  Y. Li, M. Swertz, G. Vera, J. Fu, R. Breitling, and R.C. Jansen. designGG: 
  An R-package and Web tool for the optimal design of genetical genomics 
  experiments. (submitted) \cr
  http://gbic.biol.rug.nl/designGG \cr
  Y. Li, R. Breitling and R.C. Jansen. Generalizing genetical 
     genomics: the added value from environmental perturbation, Trends Genet
     (2008) 24:518-524. \cr
  E. Wit and J. McClure. Statistics for Microarrays: Design, Analysis
      and Inference. (2004) Chichester: Wiley.    
}
\author{  
 Yang Li <yang.li@rug.nl>, Gonzalo Vera <gonzalo.vera.rodriguez@gmail.com> \cr
  Rainer Breitling <r.breitling@rug.nl>, Ritsert Jansen <r.c.jansen@rug.nl> 
  }
\seealso{ 
  \code{\link{initialDesign}}, \code{\link{designScore}}, 
  \code{\link{updateDesign}}, \code{\link{acceptanceProbability}}, \cr
  \code{\link{experimentDesignTable}}, \code{\link{plotAllScores}}
}
\examples{
  library(designGG)
  #load genotype
  data(genotype)
  #Example:  single-channel experiment with 2 environmental factors,
  #each with 5 level, and there will be four samples per condition(nTuple=4).
  optimalDesign <- designGG ( genotype, nSlides=NULL, nTuple=4, nEnvFactors=2,
                        nLevels=c(5,5),Level=NULL,  bTwoColorArray=FALSE,
                        initial=NULL, weight=1, region=NULL, optimality="A", 
                        method="SA", nIterations=100, n.search=2, endTemp=1e-10,
                        startTemp=1, maxTempStep=0.9, plotScores=TRUE, 
                        directory=NULL, fileName=NULL, envFactorNames=NULL, 
                        writingProcess=FALSE )
  #Example 2:  dual-channel experiment with 2 environmental factors,
  #each with 5 level. There are 50 slides available.
  optimalDesign <- designGG ( genotype, nSlides=50, nTuple=NULL, nEnvFactors=2,
                        nLevels=c(5,5),Level=NULL,  bTwoColorArray=TRUE,
                        initial=NULL, weight=1, region=NULL, optimality="A", 
                        method="SA", nIterations=100, n.search=2, endTemp=1e-10,
                        startTemp=1, maxTempStep=0.9, plotScores=TRUE, 
                        directory=NULL, fileName=NULL, envFactorNames=NULL, 
                        writingProcess=FALSE )
  #result
  optimalDesign$arrayDesign
  optimalDesign$conditionDesign
  plotAllScores(optimalDesign$plot.obj)
}
\keyword{ method }
