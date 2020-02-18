## Overview
The POPScomp program computes the Solvent Accessible Surface Area (SASA)
of a given PDB structure. If the structure is a complex, *i.e.* composed of more than one chain
containing protein or RNA/DNA, POPScomp creates internally all pair combinations
of chains to compute the buried SASA upon complexation. Details of those routines
are explained in the [published papers](https://github.com/Fraternalilab/POPScomp/#References)
on implicit solvent, POPS and POPSCOMP. The POPScomp program unfies the latter two methods.

There are several options to run POPScomp, refer to the [Wiki](https://github.com/Fraternalilab/POPScomp/wiki) for more details:
1. Use our POPScomp *Shiny* server at http://popscomp.org:3838 .
The interface allows for easy processing of PDB structures and shows result tabs for atom, residue, chain and molecule SASAs.
2. Download the POPScomp [Docker Image](https://github.com/Fraternalilab/POPScomp/#DockerImage) and use the *Shiny* App on your local computer without any further installation.
3. Clone the [POPScomp GitHub repository](https://github.com/Fraternalilab/POPScomp), compile POPSC and run the *Shiny* App on your local computer.


## Software
Please use the [Issues](https://github.com/Fraternalilab/POPScomp/issues) page for bug reports and add a Star to the repository to support the software maintainers.

The [POPScomp GitHub repository](https://github.com/Fraternalilab/POPScomp) contains the following components:
### POPSC  ![C CI](https://github.com/Fraternalilab/POPScomp/workflows/C%20CI/badge.svg) 
A GNU Autotools package of the POPS C program that computes SASA for a given structure.
  
### POPSR  ![R CI](https://github.com/Fraternalilab/POPScomp/workflows/R%20CI/badge.svg)
An R package to
- split complexes into single and pair components to compute buried SASA using POPSC and 
- provide a *Shiny* App as interface to the POPSC program.

### DockerImage
A *Docker* image of POPScomp has been developed and is currently (mid Jan 2020) being tested.
Watch this space for further updates.


### FunPDBe
Scripts to run POPScomp over the PDB database and feed the output into the FunPDBe project.


## Servers

### POPScomp server
* [POPScomp](http://popscomp.org:3838)

Results will be stored on the server for maximally one day.
For permanent storage, please download your results via the 'Download' buttons.

### FunPDBe
POPScomp is part of the [FunPDBe resources](https://www.ebi.ac.uk/pdbe/funpdbe/deposition).

### Note on legacy sites
* mathbio.nimr.mrc.ac.uk/POPSCOMP : This was the original address of the POPS and POPSCOMP servers in the Division of Mathematical Biology at the NIMR institute, which has been closed in 2015 and became a part of the Francis Crick Institute.
* mathbio.crick.mrc.ac.uk/POPSCOMP : For a short while, the POPS and POPSCOMP servers were running in the Francis Crick Institute, but the institute's administration decided to stop support.


## References
Users publishing results obtained with the program and its applications
should acknowledge its use by citation.

### Implicit solvent
Fraternali, F. and van Gunsteren, W.F.<br>
*An efficient mean solvation force model for use in molecular dynamics simulations of proteins in aqueous solution.*<br>
**Journal of Molecular Biology** 256 (1996) 939-948.<br>
[![doi](https://img.shields.io/badge/doi-10.1016%2Fj.jmb.2014.03.010-blue.svg?style=flat)](https://dx.doi.org/10.1016%2Fj.sbi.2014.04.003)
[![pubmed](https://img.shields.io/badge/pubmed-24681267-blue.svg?style=flat)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4045398/)

### POPS method
Fraternali, F. and Cavallo, L.<br>
*Parameter optimized surfaces (POPS): analysis of key interactions and conformational changes in the ribosome.*<br>
**Nucleic Acids Research** 30 (2002) 2950-2960.<br>
[![doi](https://img.shields.io/badge/doi-10.1016%2Fj.jmb.2014.03.010-blue.svg?style=flat)](https://dx.doi.org/10.1093%2Fnar%2Fgkf373)
[![pubmed](https://img.shields.io/badge/pubmed-24681267-blue.svg?style=flat)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC117037/)

### POPS server
Cavallo, L., Kleinjung, J. and Fraternali, F.<br>
*POPS: A fast algorithm for solvent accessible surface areas at atomic and residue level.*<br>
**Nucleic Acids Research** 31 (2003) 3364-3366.<br>
[![doi](https://img.shields.io/badge/doi-10.1016%2Fj.jmb.2014.03.010-blue.svg?style=flat)](https://dx.doi.org/10.1093%2Fnar%2Fgkg601)
[![pubmed](https://img.shields.io/badge/pubmed-24681267-blue.svg?style=flat)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC169007/)

### POPSCOMP server
Kleinjung, J. and Fraternali, F.<br>
*POPSCOMP: an automated interaction analysis of biomolecular complexes.*<br>
**Nucleic Acids Research** 33 (2005) W342-W346.<br>
[![doi](https://img.shields.io/badge/doi-10.1016%2Fj.jmb.2014.03.010-blue.svg?style=flat)](https://dx.doi.org/10.1093%2Fnar%2Fgki369)
[![pubmed](https://img.shields.io/badge/pubmed-24681267-blue.svg?style=flat)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1160130/)


## License and Copyright
Usage of the software and server is free under the GNU General Public License v3.0.

### Copyright Holders, Authors and Maintainers 
- 2002-2019 Franca Fraternali (author, maintainer) franca.fraternali@kcl.ac.uk
- 2008-2019 Jens Kleinjung (author, maintainer) jens@jkleinj.eu

### Contributors
- 2002 Kuang Lin and Valerie Hindie (translation to C)
- 2002 Luigi Cavallo (parametrisation)

