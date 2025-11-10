## POPScomp: Solvent Accessible Surface Areas of Biomolecules and their Complexes
[![release](https://img.shields.io/badge/release-v3.2.1-green?logo=github)](https://github.com/Fraternalilab/POPScomp/releases)
[![release](https://img.shields.io/badge/release-v3.2.1-green?logo=docker)](https://hub.docker.com/u/popscomp)
[![wiki](https://img.shields.io/badge/wiki-orange)](https://github.com/Fraternalilab/POPScomp/wiki)
[![License](https://img.shields.io/badge/license-GPL-blue.svg)](http://www.gnu.org/licenses/gpl-3.0.en.html)
[![DOI](https://zenodo.org/badge/182454886.svg)](https://zenodo.org/badge/latestdoi/18245488?event=push6)
[![C CI](https://github.com/Fraternalilab/POPScomp/workflows/C%20CI/badge.svg)](https://github.com/thereapr/reincarnated_mod/actions/workflows/ccpp.yml)
[![R CI](https://github.com/Fraternalilab/POPScomp/workflows/R%20CI/badge.svg)](https://github.com/thereapr/reincarnated_mod/actions/workflows/main.yml)


The POPScomp program computes the Solvent Accessible Surface Area (SASA)
of a given PDB structure. If the structure is a complex, *i.e.* composed of more than one chain
containing protein or RNA/DNA, POPScomp computes additionally the SASA buried between chains.
The name 'POPS' is derived from **P**arameter **OP**timised **S**urfaces, because its parametrisation
was performed by regression against a large set of protein and RNA/DNA structures.
Details of those routines are explained in the
[published papers](https://github.com/Fraternalilab/POPScomp/#References).
(Note: POPScomp embodies both methods, POPS and POPSCOMP).

## Installation
There are 2 ways to (install and) use POPScomp,
please refer to the [Wiki](https://github.com/Fraternalilab/POPScomp/wiki) for details:
1. Download the [POPScomp *Docker* image](https://hub.docker.com/repository/docker/jkleinj150999/popscomp/)
  (for backward compatibility, this is the [legacy POPScomp *Docker* image](https://hub.docker.com/u/popscomp))
  and run the *Shiny* App on your local computer.
3. To install from source code, either download the [latest release](https://github.com/Fraternalilab/POPScomp/releases/latest) or the development version by cloning this POPScomp repository.
[Compile/install](https://github.com/Fraternalilab/POPScomp/wiki/Installation) the program suite and run the *Shiny* App on your local computer.


## Usage
Enter a PDB identifier or upload a custom PDB file and press 'run POPScomp'.
Results will be shown as tables for different resolutions: Atom, Residue, Chain and Molecule.
Download individual results *via* the buttons underneath the tables or bundled *via* the
'Download All Results' button.

![POPScomp Shiny Interface](https://github.com/Fraternalilab/POPScomp/blob/master/POPSR/inst/popsr/png/POPScomp_Shiny.png)


## Documentation
Please refer to the [Wiki](https://github.com/Fraternalilab/POPScomp/wiki) for documentation.


## Software Design
This POPScomp repository contains the following components:

#### POPSC
A GNU Autotools package of the POPS C program that computes SASA for a given structure.
  
#### POPSR
An R package to
- split complexes into single and pair components to compute buried SASA using POPSC and 
- provide a *Shiny* App as interface to the POPSC program.

#### *Docker* Image
A [POPScomp *Docker* image](https://hub.docker.com/repository/docker/jkleinj150999/popscomp/) is available that provides 
the program suite in fully functional form,
*i.e.* pre-compiled POPSC, RStudio with POPSR and all dependencies installed.
Please refer to the [Wiki](https://github.com/Fraternalilab/POPScomp/wiki/Docker-Image)
for detailed guidelines on downloading and usage of the Docker image.

#### FunPDBe
A *Make* pipeline to run POPScomp over the PDB database and to feed the JSON output into the [Elixir FunPDBe project](https://www.ebi.ac.uk/pdbe/pdbe-kb/funpdbe).
To include newly added PDB structures, the pipeline is run *via* a cron *@weekly* job on a virtual
Ubuntu machine hosted on a Synology DS218+ computer. 
This directory has been included purely as a reference for the FunPDBe project and is not intended as a user-facing application.
However, some users might find the pipeline or parts of it useful for their own purposes.


## Servers

#### Legacy POPScomp *Shiny* server at [popscomp.org](http://popscomp.org:3838)
This server has been discontinued in November 2025.
Please download the Docker image instead to run POPScomp on your local machine.

#### FunPDBe
POPScomp is part of the [FunPDBe resources](https://www.ebi.ac.uk/pdbe/funpdbe/deposition).

#### Note on older legacy sites
* mathbio.nimr.mrc.ac.uk/POPSCOMP : This was the original address of the POPS and POPSCOMP servers in the Division of Mathematical Biology at the NIMR institute, which has been closed in 2015 and became a part of the Francis Crick Institute.
* mathbio.crick.mrc.ac.uk/POPSCOMP : Only for a short while the POPS and POPSCOMP servers were supported by the Francis Crick Institute.


## References
Users publishing results obtained with the program and its applications
should acknowledge its use by citation.

#### POPScomp software
[![DOI](https://zenodo.org/badge/182454886.svg)](https://zenodo.org/badge/latestdoi/182454886)

#### Implicit solvent
Kleinjung, J. and Fraternali, F.<br>
*Design and Application of Implicit Solvent Models in Biomolecular Simulations.*<br>
**Current Opinion in Structural Biology** 25 (2014) 126-134.<br> 
[![doi](https://img.shields.io/badge/doi-10.1016/j.sbi.2014.04.003-blue.svg?style=flat)](https://dx.doi.org/10.1016/j.sbi.2014.04.003)
[![pubmed](https://img.shields.io/badge/pubmed-4045398-blue.svg?style=flat)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4045398/)

Fraternali, F. and van Gunsteren, W.F.<br>
*An efficient mean solvation force model for use in molecular dynamics simulations of proteins in aqueous solution.*<br>
**Journal of Molecular Biology** 256 (1996) 939-948.<br>
[![doi](https://img.shields.io/badge/doi-10.1016%2Fj.jmb.2014.03.010-blue.svg?style=flat)](https://dx.doi.org/10.1016%2Fj.sbi.2014.04.003)
[![pubmed](https://img.shields.io/badge/pubmed-24681267-blue.svg?style=flat)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4045398/)

#### POPS method
Fraternali, F. and Cavallo, L.<br>
*Parameter optimized surfaces (POPS): analysis of key interactions and conformational changes in the ribosome.*<br>
**Nucleic Acids Research** 30 (2002) 2950-2960.<br>
[![doi](https://im)](https://dx.doi.org/10.1093%2Fnar%2Fgkf373)
[![pubmed](https://img.shields.io/badge/pubmed-117037-blue.svg?style=flat)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC117037/)

#### POPS server
Cavallo, L., Kleinjung, J. and Fraternali, F.<br>
*POPS: A fast algorithm for solvent accessible surface areas at atomic and residue level.*<br>
**Nucleic Acids Research** 31 (2003) 3364-3366.<br>
[![doi](https://img.shields.io/badge/doi-0.1093%2Fnar%2Fgkg601-blue.svg?style=flat)](https://dx.doi.org/10.1093%2Fnar%2Fgkg601)
[![pubmed](https://img.shields.io/badge/pubmed-169007-blue.svg?style=flat)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC169007/)

#### POPSCOMP method and server
Kleinjung, J. and Fraternali, F.<br>
*POPSCOMP: an automated interaction analysis of biomolecular complexes.*<br>
**Nucleic Acids Research** 33 (2005) W342-W346.<br>
[![doi](https://img.shields.io/badge/doi-10.1093%2Fnar%2Fgki369-blue.svg?style=flat)](https://dx.doi.org/10.1093%2Fnar%2Fgki369)
[![pubmed](https://img.shields.io/badge/pubmed-1160130-blue.svg?style=flat)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1160130/)


## Contributing
Please feel free to open an issue on the [Issues](https://github.com/Fraternalilab/POPScomp/issues) page if you have any difficulty with POPScomp.
We have no formal procedure to process contributed code, but we will revise code submitted *via* pull requests and will acknowledge authors accordingly.


## License and Copyright
Usage of the software and server is free under the GNU General Public License v3.0.

#### Copyright Holders, Authors and Maintainers 
- 2002-2021 Franca Fraternali (author, maintainer) franca.fraternali@kcl.ac.uk
- 2008-2021 Jens Kleinjung (author, maintainer) jens@jkleinj.eu

#### Contributors
- 2002 Kuang Lin and Valerie Hindie (translation to C)
- 2002 Luigi Cavallo (parametrisation)

