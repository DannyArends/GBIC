Metanetwork: Genetic study of metabolites and network reconstruction
====================================================================
MetaNetwork maps metabolite quantitative trait loci (mQTLs) underlying variation in metabolite 
abundance in individuals of a segregating population using a two-part model to account for the 
nature of metabolite data.

This model combines the analysis of the binary traits (positive/not-available) with conditional 
analysis of the quantitative trait (numeric) among individuals with a positive binary phenotype. 
Simulation procedures are used to assess statistical significance 

MetaNetwork predicts the network of potential associations between metabolites using correlations 
of mQTL profiles or abundance profiles 

MetaNetwork generates files of predicted networks, which can be visualized using Cytoscape, and 
optionally relates multiple mass peaks per metabolite that may be consequence of isotopes or 
charge difference.

Analysis of about 24 metabolites takes a few minutes on a desktop computer (Pentium 4). Analysis 
of a metabolome of about 2000 metabolites will take around four days.

Dependencies
------------
R 2.10.0 or higher with the R package 'qvalue' installed.

Installation
------------
Download the provided .zip R package (Be carefull not to download the source) for windows 32/64 bit, or 
the tar.gz for other operating systems.
On the commandline use: R CMD INSTALL MetaNetwork_X.X-X.zip or R CMD INSTALL MetaNetwork_X.X-X.tar.gz

References
----------
Fu J, Swertz MA, Keurentjes JJB, Jansen RC. MetaNetwork: a computational tool for the genetic study of metabolism. Nature Protocols (2007). 

Disclaimer
----------
This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License,
version 3, as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but without any warranty; without even the implied warranty of
merchantability or fitness for a particular purpose.  See the GNU
General Public License, version 3, for more details.

A copy of the GNU General Public License, version 3, is available
at [http://www.r-project.org/Licenses/GPL-3](http://www.r-project.org/Licenses/GPL-3 "GPL-3 Licence")
Copyright (c) 2007-2011 
Jingyuan Fu <j.fu@rug.nl>, 
Morris Swertz <m.a.swertz@rug.nl>, 
Ritsert Jansen <r.c.jansen@rug.nl> 

