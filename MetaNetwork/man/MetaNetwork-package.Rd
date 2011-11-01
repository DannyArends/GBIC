\name{ Overview }
\alias{ MetaNetwork-package }
\title{ MetaNetwork package }
\description{
  MetaNetwork maps metabolite quantitative trait loci (mQTLs) underlying 
  variation in metabolite abundance in individuals of a segregating population 
  with a two-part model to account for the often observed spike in the 
  distribution of metabolite data. MetaNetwork predicts and visualizes 
  potential associations between metabolites using correlations of mQTL 
  profiles, rather than of abundance profiles. Simulation and permutation 
  procedures are used to assess statistical significance.
}
\note{ 
  To install MetaNetwork package:
  Download the MetaNetwork package from: \url{http://gbic.biol.rug.nl/supplementary/2007/MetaNetwork}
  (Windows/RGui) click 'Packages'; 'install package(s) from local zip files'; browse to the MetaNetwork zip file.
  
  MetaNetwork requires the \code{qvalue} package for false discovery rate control:
  (Windows/RGui) click 'Packages'; 'Install package(s)...'; select CRAN mirror; select qvalue.
  
  The developers of MetaNetwork thank Gonzalo Vera, Bruno Tesson and Richard
  Scheltema for helping to develop this R-package, and Martijn Dijkstra and 
  Rainer Breitling for helping to predict multiple peaks belonging to the same 
  metabolite
}
\references{ 
  Please cite publication below when reporting results produced using this package: 
  Fu J, Swertz MA, Keurentjes JJB, Jansen RC. MetaNetwork: a computational tool for the 
  genetic study of metabolism. Nature Protocols (2007). 
}
\author{ 
Jingyuan Fu <j.fu@rug.nl>, 
Morris Swertz <m.a.swertz@rug.nl>,
Joost Keurentjes <Joost.Keurentjes@wur.nl>, 
Ritsert Jansen <r.c.jansen@rug.nl> }