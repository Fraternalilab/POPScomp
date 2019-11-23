## Overview
The POPS program computes the Solvent Accessible Surface Area (SASA)
of a given PDB structure. If the structure is composed of more than one chain
containing protein or RNA/DNA, POPScomp creates internally all pair combinations
of chains to compute the buried SASA upon complexation. Details of those routines
are explained in the published papers on POPS and POPSCOMP.

POPScomp (*Shiny* app) shows tabs for atom, residue, chain and molecule SASAs.
The tables are initialised without any values and therefore the user sees
the table header and below the notice 'No data available in table'.
After selecting a PDB identifier or file and pressing 'run POPScomp',
the Shiny sever runs the POPS program on components of the PDB file
and the tables automatically refresh to show the resulting SASA values.
Because running POPS is a system call, the success of the computation
is returned as exit code and shown below the 'run POPScomp' button:
* 0 - Success
* 1 - Catchall for general errors
* 2 - Misuse of shell builtins (according to Bash documentation)
* 126 - Command invoked cannot execute
* 127 - Command not found
* 128 - Invalid argument to exit
* 128+n - Fatal error signal 'n'
* 130 - Script terminated by Control-C
* 255\* - Exit status out of range


## Software

### Packages
Since POPScomp 3.0 (04.2019), the modules *POPSC*, *POPSR*, *FunPDBe* and
*Docker* are part of the POPScomp project.
1. *POPSC*: A GNU Autotools package of the POPS C program that computes SASA
  for a given structure.
2. *POPSR*: An R package containing an R program that
    - splits complexes into single and pair components to compute buried SASA
	  using POPSC and 
    - provides a Shiny server as interface to the POPScomp (POPSC and POPSR) programs.
3. *FunPDBe*: Scripts to run POPScomp over the PDB database and feed the output into
  the FunPDBe project.
4. *Docker*: A Docker image is under development and scheduled for version 3.1.

### Source code
Fraternali lab's
[POPScomp GitHub page](https://github.com/Fraternalilab/POPScomp).
Please use that site for bug reports and add a Star to the repository
to support the software maintainers.

### Docker image
A *Docker* image of POPScomp can be pulled from the Google Cloud.
```
docker pull eu.gcr.io/high-hue-217311/popscomp:0.1
```

## Servers

### POPScomp server
* [POPScomp](http://popscomp.org:3838)

Results will be stored on the server for maximally one day.
For permanent storage, please download your results via the 'Download' buttons.

### FunPDBe
POPScomp is part of the [FunPDBe resources](https://www.ebi.ac.uk/pdbe/funpdbe/deposition).


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
Usage of the software and server is free, the code license is GPL3.

### Copyright Holders, Authors and Maintainers 
- 2002-2019 Franca Fraternali (author, mainatainer) franca.fraternali@kcl.ac.uk
- 2008-2019 Jens Kleinjung (author maintainer) jens@jkleinj.eu

### Contributors
- 2002 Kuang Lin and Valerie Hindie (translation to C)
- 2002 Luigi Cavallo (parametrisation)

