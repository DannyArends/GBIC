AffyGG: computational protocols for genetical genomics with Affymetrix arrays
====================================================================
Affymetrix arrays use multiple probes per gene to measure mRNA 
abundances. Standard software takes averages over probes. Important 
information may be lost if polymorphisms in the mRNA affect the 
hybridization of individual probes.

AffyGG is software to analyze genetical genomics experiments in human, mouse 
and other organisms: (i) an R package providing functions for QTL analysis 
at the individual probe level and (ii) Perl scripts providing custom tracks 
in the UCSC Genome Browser to check for sequence polymorphisms in probe regions.

Dependencies
------------
R 2.10.0 or higher with the R package 'affy' and 'nlme' installed.

Installation
------------
Download the provided .zip R package (Be carefull not to download the source) for windows 32/64 bit, or 
the tar.gz for other operating systems.
On the commandline use: R CMD INSTALL affyGG_X.X-X.zip or R CMD INSTALL affyGG_X.X-X.tar.gz

References
----------
Alberts R, Vera G, Jansen RC. AffyGG: computational protocols for genetical genomics with Affymetrix arrays - Bioinformatics (2008).

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
Copyright (c) 2008-2011 
Rudi Alberts
Gonzalo Vera
Ritsert Jansen <r.c.jansen@rug.nl> 

