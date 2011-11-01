\name{MetaNetwork}
\alias{MetaNetwork}
\title{genetic study of metabolites and network reconstruction }
\description{
 An integrated computational protocol to run a complete genetic analysis workflow on 
 metabolites using diverse MetaNetwork methods for quantitative trait analysis,
 network reconstruction and Cytoscape network visualization. 
}
\usage{
MetaNetwork( markers, genotypes, traits, spike, qtlProfiles = NULL, 
    peaks = NULL, qtlThres = NULL, qtlSumm = NULL, corrZeroOrder = NULL, 
    corrSecondOrder = NULL, corrMethod = "qtl", corrThres = 0, 
    cytoFiles = TRUE, outputdir = "./MetaNetwork")
}
\details{
  First, MetaNetwork maps metabolite quantitative trait loci (mQTLs) underlying 
  variation in metabolite abundance in individuals of a segregating population 
  using a two-part model to account for the nature of 
  metabolite data (step A, \code{\link{qtlMapTwoPart}}). This model combines the analysis of the binary traits 
  (positive/not-available) with conditional analysis of the quantitative 
  trait (numeric) among individuals with a positive binary phenotype. Simulation
  procedures are used to assess statistical significance (step B, \code{\link{qtlThreshold}}, 
  \code{\link{qtlFDR}}). MetaNetwork will summarize the information about significant
  mQTLs (step C, \code{\link{qtlSummary}}).

  Then, MetaNetwork predicts the network of potential associations 
  between metabolites using correlations of mQTL profiles or 
  abundance profiles (step D, \code{\link{qtlCorrZeroOrder}}; 
  step E, \code{\link{qtlCorrSecondOrder}}). 
  Optionally, permutation procedures can be used to assess statistical 
  significance (step F, \code{\link{qtlCorrThreshold}}). 
  
  Finally, MetaNetwork generates files of predicted networks, which can 
  be visualized using Cytoscape (step G, \code{\link{createCytoFiles}}), and
  optionally relates multiple mass peaks per metabolite that may be consequence
  of isotopes or charge difference (step H, \code{\link{findPeakMultiplicity}}). 
  
  Analysis of about 24 metabolites takes a few minutes on a 
  desktop computer (Pentium 4). Analysis of a metabolome of about 2000 
  metabolites will take around four days. In addition, MetaNetwork 
  is able to integrate high-throughput data from future metabolomics, 
  transcriptomics and proteomics experiments in conjunction with phenotypic 
  data.
  
  After running MetaNetwork with defaults, the R console will show:
  \preformatted{
>MetaNetwork (markers=markers, genotypes=genotypes, traits=traits, spike=4, qtlThres=3.79)        
Step A: QTL mapping.... 
         result in R object 'qtlProfiles' 
         result in ./MetaNetwork/qtlProfiles.csv 
         process time 29.25 sec 

Step B: Simulation test for QTL significance threshold....skipped 
         using user-provided QTL threshold: 3.79 

Step C: QTL summary.... 
         result in R object: 'qtlSumm' 
         result in ./MetaNetwork/qtlSumm.csv 
         process time 1.66 sec 

Step D: Zero-order correlation .... 
         result in R object: 'corrZeroOrder' 
         result in ./MetaNetwork/corrZeroOrder.csv 
         process time 2.97 sec 

Step E: 2nd-order correlation .... 
         result in R object: 'corrSecondOrder' 
         result in ./MetaNetwork/corrSecondOrder.csv 
         process time 9.58 sec 

Step F: Permutation test for 2nd-order correlation significance threshold...skipped 
         using user-provided correlation threshold: 0 

Step G: Create Cytoscape network files... 
         SIF file is: ./MetaNetwork/network.sif 
         EDA file is: ./MetaNetwork/network.eda 

Step H: Find peak multiplicity...skipped
  }  
}
\arguments{
  \item{markers}{ matrix of markers (rownames) and their chromosome numbers 
  (column 1) and centi-Morgan positions (cM, column 2), ordered by position.
  See \code{\link{markers}} example data.}
  \item{genotypes}{ matrix of genotypes for each marker (rownames) and individual 
  (columnnames), as numeric values 1, 2 or NA when missing.
  See \code{\link{genotypes}} example data.}
  \item{traits}{ matrix of phenotypes for each trait (rownames) and individual 
  (columnnames), as numeric or NA when missing.
  See \code{\link{traits}} or \code{\link{traits2}} example data.}
  \item{spike}{ numeric cut-off value to separate absent (qualitative) from 
  available (quantitative) trait abundance.}
  \item{qtlProfiles}{ (optional) matrix of QTL mapping of traits (rownames) 
  to markers (columnnames), as \eqn{-log_{10}(p)} values. 
  If \code{qtlProfiles} is NULL (default), MetaNetwork will call function 
  \code{\link{qtlMapTwoPart}} to generate the \code{qtlProfiles}. 
  Otherwise, MetaNetwork will use the provided QTL results. 
  See \code{\link{qtlProfiles}} example data.}
  \item{qtlThres}{ (optional) numeric \eqn{-log_{10}(p)} threshold value for significant QTLs.  
  If \code{qtlThres} is NULL (default), the QTL significance threshold will be estimated 
  by simulation using \code{\link{qtlThreshold}} at \code{alpha = 0.05} and 
  \code{n.simulations = 1000}. Furthermore, the QTL significance threshold is also estimated by
  controlling the false discovery rate  \code{fdrThres = 0.05} 
  using \code{\link{qtlFDR}}. 
  By default, the most stringent outcome of either \code{\link{qtlThreshold}} or \code{\link{qtlFDR}} is used. 
  Otherwise, MetaNetwork will use provided threshold. }
  \item{qtlSumm}{ (optional) data.frame with the summary of each QTL.
  If \code{qtlSumm} is NULL (default), MetaNetwork will 
  call function \code{\link{qtlSummary}} to summarize QTL effects.
  See \code{\link{qtlSumm}} example data.}
  \item{corrZeroOrder}{ (optional) the matrix of zero-order correlation coefficient between 
  metabolites.  
  If \code{corrZeroOrder} is NULL (default), zero-order correlation coefficient 
  will be calculated for QTL profiles using function \code{\link{qtlCorrZeroOrder}}.
  See \code{\link{corrZeroOrder}} example data.}
  \item{corrSecondOrder}{ (optional) matrix of second-order partial correlation between 
  metabolites.
  If \code{corrSecondOrder} is NULL (default), second-order partial correlation 
  will be calculated for QTL profiles using function \code{\link{qtlCorrSecondOrder}}.
  See \code{\link{corrSecondOrder}} example data.}
  \item{corrMethod}{ (optional) character string indicating which correlation method, 
  either "qtl" or "abundance".  
  If \code{corrMethod} is "qtl" (default), MetaNetwork 
  will call function \code{\link{qtlCorrZeroOrder}} to calculate the correlation between 
  QTL profiles. Otherwise, when \code{corrMethod} is "abundance", MetaNetwork will 
  use Spearman corrrelation via function \code{\link{cor}} to calculate the correlation between metabolite 
  abundance profiles.}
  \item{corrThres}{ (optional) numeric threshold for significant partial correlation coefficients.  
  If \code{corrThres} is NULL, the empirical threshold is estimated by permutation 
  using function \code{qtlCorrThreshold} with \code{n.permutations} = 10000. 
  Otherwise, the provided threshold is used. Default is 0. }
  \item{peaks}{ (optional) matrix of mass/charge peaks (column1) for each trait (rownames).
  If \code{peaks} is set, MetaNetwork will call \code{\link{findPeakMultiplicity}} to
  relate multiple mass peaks for correlated traits. 
  See \code{\link{peaks2}} example peaks data for unidentified metabolite example traits data \code{\link{traits2}}.}
  \item{cytoFiles}{ (optional) boolean value that indicates if files for network visualization 
  in Cytoscape should be created. 
  If TRUE (default) MetaNetwork will call function \code{\link{createCytoFiles}} 
  to create two network files in \code{outputdir} for the significant correlations amongst metabolites: 
  'network.sif' and 'network.eda'.  
  }
  \item{outputdir}{ (optional) output directory where generated data files will be stored. 
  Default is "./MetaNetwork" }
}
\value{
  \item{qtlProfiles}{matrix of QTL mapping of traits (rownames) 
  to markers (columnnames) as log-transformed "p values" [ \eqn{-log_{10}(p)}],  
  see \code{\link{qtlMapTwoPart}}. A +/- sign 
  is added to indicate the direction of the additive effect:                    
  values are positive if the QTL has higher metabolite abundance for individuals       
  carrying the genotype 2 than those carrying the genotype 1; values are negative 
  otherwise. 
  See \code{\link{qtlProfiles}} example data.}                                                     
  \item{qtlThres }{estimated QTL significance threshold.
  See function \code{\link{qtlThreshold}}.}
  \item{qtlSumm }{ data frame with QTL summary. 
  See \code{\link{qtlSumm}} example data. } 
  \item{corrZeroOrder }{matrix of zero order correlation of QTL profiles.
  See \code{\link{corrZeroOrder}} example data.}
  \item{corrSecondOrder }{matrix of 2nd order correlation of QTL profiles.
  See \code{\link{corrSecondOrder}} example data.}
  \item{corrPermutations }{vector of the permutations of maximum, absolute correlation values.
  See function \code{\link{qtlCorrThreshold}}.}  
  \item{corrThres }{numeric correlation threshold.
  See function \code{\link{qtlCorrThreshold}}.}
  \item{cytoFiles}{network files "network.sif" and "network.eda" for cytoscape are produced in \code{outputdir}. 
  See function \code{\link{createCytoFiles}}.} 
  \item{peakMultiplicity}{If \code{peaks} is not NULL, data frame with Peak Multiplicity summary.
  See \code{\link{peakMultiplicity}} example data.}
  \item{resultFiles}{If \code{outputdir} is not NULL, the above outputs will be also saved in 
  files "qtlProfiles.csv", "qtlSumm.csv", "corrZeroOrder.csv","corrSecondOrder.csv", 
  "corrPermutations.csv", "peakMultiplicity.csv", respectively. A summary of 
  analysis processing, results objects and output files can be seen in the R console
  and is saved in file "output.txt".}  
}
\references{ 
  Fu J, Swertz MA, Keurentjes JJB, Jansen RC. MetaNetwork: a computational tool for the 
  genetic study of metabolism. Nature Protocols (2007). 
  
  \url{http://gbic.biol.rug.nl/supplementary/2007/MetaNetwork}
}
\author{ 
Jingyuan Fu <j.fu@rug.nl>, 
Morris Swertz <m.a.swertz@rug.nl>, 
Ritsert Jansen <r.c.jansen@rug.nl>
}
\note{ 
The names of individuals (columnnames) must be consistent over \code{genotypes} and \code{traits}. 
The names of peaks (rownames) must be consistent over \code{peaks} and \code{traits}.
}
\seealso{ 
  Use \code{\link{markers}}, \code{\link{genotypes}} and \code{\link{traits}} as example data sets or use \code{\link{loadData}} to load your own data.
  Use \code{\link{qtlMapTwoPart}} for the calculation of \code{qtlProfiles}.
  Use \code{\link{qtlThreshold}} and \code{\link{qtlFDR}} for the estimation of \code{qtlThres} QTL significance threshold.
  Use \code{\link{qtlCorrZeroOrder}} and \code{\link{qtlCorrSecondOrder}} for the calculation of zero order and second order correlation for \code{corrZeroOrder} and \code{corrSecondOrder} respectively.
  Use \code{\link{qtlCorrThreshold}} for the estimation of \code{corrThres} correlation significance threshold.
  Use \code{\link{qtlSummary}} for the generation of \code{qtlSumm} QTL summary. 
  Use \code{\link{createCytoFiles}} for the generation of Cytoscape network files. 
  Use \code{\link{findPeakMultiplicity}} for the relation of isotopic or differentially charged metabolites.
}
\examples{
## load the example data provided with this package
data(genotypes)
data(traits)
data(markers)

#set qtlThres
qtlThres    <- 3.79

#run metanetwork with predefined thresholds
\dontrun{
MetaNetwork (markers=markers, genotypes=genotypes, traits=traits, spike=4, qtlThres=qtlThres)

##OR: load data from csv
#genotypes <- loadData("genotypes.csv")
#traits    <- loadData("traits.csv")
#markers   <- loadData("markers.csv")
#MetaNetwork (markers=markers, genotypes=genotypes, traits=traits2, 
#             qtlThres=qtlThres, spike=4) 
             
##OR: let MetaNetwork estimate qtlThres and identify multiple peaks
#data(genotypes)
#data(traits2)
#data(markers)
#data(peaks2)
#MetaNetwork (markers=markers, genotypes=genotypes, traits=traits2, 
#             peaks=peaks2, spike=4)                           
  
##show part of the qtlProfiles
qtlProfiles[1:5,1:5]

##show part of the qtl summary
qtlSumm[1:5,]

##show part of the zero order correlation
corrZeroOrder[1:5,1:5]

##show part of the second order correlation
corrSecondOrder[1:5,1:5]

##plot the qtlProfiles
qtlPlot(markers, qtlProfiles, qtlThres)  
}
##load network.sif and network.eda into Cytoscape 
}