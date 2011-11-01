\name{plotAllScores} 
\alias{plotAllScores} 
\title{Plot scores profiles}

\description{
   plot all scores and cooling at each iteration during simulated annealing 
   process
} 
\usage{ plotAllScores(plot.obj,fileName=NULL) 
} 
\arguments{     
  \item{plot.obj}{ An list containing: scores, cooling, startTemp, temperature, 
                       temperature.step, nIterations and optimality. 
                       Details can be found below. 
  }
  \item{scores}{
	 A- or D- optimality score of all accepted designs during optimization
	   process.
	}
  \item{cooling}{ 
	descibes the colling step in the Simulated Annealing process
	cooling= (new.score - now.score)/ now.score
	}
 \item{startTemp}{
        starting temperature of simulated annealing process.
  }
 \item{temperature}{
	final temperature that the Simulated annealing reaches.
	}

  \item{temperatureStep}{
        temperature decreasing step for simulated annealing process.
  }		
 \item{nIterations}{
        Number of iterations of the simulated annealing method.
  }
 \item{optimality}{
   type of optimality, i.e. "A" (A-optimality) or "D"
      (D-optimality). A-optimality minimizes $Trace((X'X)^{-1})$,
      which corresponds to minimum average variance of the
      parameter estimates. D-optimality minimizes $det(X'X)^{-1}$,
      which corresponds to minimum generalized variance of the
      parameter estimates.
  }
 \item{fileName}{
       A plot(in png format) of the all scores during SA process will be 
       produced. The users can specify plot name by setting fileName.
       If NULL (the default) it produces files starting with "myDesignGG".
  }
}
 

\value{
  Draw a plot that visualizeds the socres (y-axis) at each iteration during 
  the simulated annealing process (x-axis is time of moving)
 
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
    The calculation of score is dependent on the choice of optimality.\cr
    Cooling is defined as (newScore-nowScore)/nowScore. 
} 
\examples{
    data(exmaplePlotObj)
    plotAllScores(plot.obj) 
}
\keyword{method}
