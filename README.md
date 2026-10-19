## POPScomp: Solvent Accessible Surface Areas of Biomolecules and their Complexes
[![release](https://img.shields.io/badge/release-v3.5.1-green?logo=github)](https://github.com/Fraternalilab/POPScomp/releases)
[![release](https://img.shields.io/badge/release-v3.5.1-green?logo=docker)](https://hub.docker.com/repository/docker/jkleinj150999/popscomp/)
[![wiki](https://img.shields.io/badge/wiki-orange)](https://github.com/Fraternalilab/POPScomp/wiki)
[![License](https://img.shields.io/badge/license-GPL-blue.svg)](http://www.gnu.org/licenses/gpl-3.0.en.html)
[![DOI](https://zenodo.org/badge/182454886.svg)](https://doi.org/10.5281/zenodo.3679109)
[![C CI](https://github.com/Fraternalilab/POPScomp/workflows/C%20CI/badge.svg)](https://github.com/Fraternalilab/POPScomp/blob/main/.github/workflows/ccpp.yml)
[![R CI](https://github.com/Fraternalilab/POPScomp/workflows/R%20CI/badge.svg)](https://github.com/Fraternalilab/POPScomp/blob/main/.github/workflows/main.yml)

The POPScomp program computes the Solvent Accessible Surface Area (SASA) of a protein or DNA/RNA molecule,
provided as PDB or CIF structure. If the structure is a complex, *i.e.* composed of more than one chain
containing protein or RNA/DNA, POPScomp computes additionally the SASA buried between chains.
The name 'POPS' is derived from **P**arameter **OP**timised **S**urfaces, because its parametrisation
was performed by regression against a large set of protein and RNA/DNA structures.
Details of those routines are explained in the
[published papers](https://github.com/Fraternalilab/POPScomp/#References).
(Note: The POPScomp release in this repository combines both methods originally published as POPS and POPSCOMP).


## Installation
There are 2 ways to (install and) use POPScomp,
please refer to the [Wiki](https://github.com/Fraternalilab/POPScomp/wiki) for details:

1. Download the latest [POPScomp *Docker* image](https://hub.docker.com/repository/docker/jkleinj150999/popscomp/tags/latest)
  and run the *Shiny* App on your local computer. The Docker image provides the program suite in fully functional form,
*i.e.* a Shiny server with integrated POPSC/POPSR as computational backbone.
Please refer to the [Wiki](https://github.com/Fraternalilab/POPScomp/wiki/Docker-Image)
for detailed guidelines on downloading and usage of the Docker image.

2. To install from source code, either download the [latest POPScomp release](https://github.com/Fraternalilab/POPScomp/releases/latest) or the development version by cloning this POPScomp repository.
[Compile/install](https://github.com/Fraternalilab/POPScomp/wiki/Installation) the program suite and run the *Shiny* App on your local computer.
```
git clone --recurse-submodules https://github.com/Fraternalilab/POPScomp.git 
```

## Usage
Enter a PDB database identifier or upload a local PDB|CIF file and press 'run POPScomp'.
SASA values and related information will be shown as tables for different resolutions: Atom, Residue, Chain and Molecule.
Download individual results *via* the buttons underneath the tables or bundled *via* the
'Download All Results' button.


## Documentation
The [Wiki](https://github.com/Fraternalilab/POPScomp/wiki) contains additional documentation about the POPS model and implementation.


## Software Design
This POPScomp repository contains the following components:

#### POPSC
A GNU Autotools C package that creates the 'pops' program, a stand-alone application for the SASA calculation of biomolecular structures.
  
#### POPSR
An R package to
- split complexes into single and pair components to compute buried SASA using 'pops' (see above) and 
- provide a *Shiny* App as interface to run POPScomp analyses and visualise the results.


## Software Sites and Releases

- [popscomp.org](http://popscomp.org) points to this GitHub repository
- [POPScomp images in Docker](https://hub.docker.com/repository/docker/jkleinj150999/popscomp)
- [POPScomp releases in GitHub](https://github.com/Fraternalilab/POPScomp/releases)
- [POPScomp releases in Zenodo](https://doi.org/10.5281/zenodo.3679109)


## References

#### Implicit solvent
Kleinjung, J. and Fraternali, F.<br>
*Design and Application of Implicit Solvent Models in Biomolecular Simulations.*<br>
**Current Opinion in Structural Biology** 25 (2014) 126-134.<br> 
[![doi](https://img.shields.io/badge/doi-10.1016%2Fj.sbi.2014.04.003-blue.svg?style=flat)](https://doi.org/10.1016/j.sbi.2014.04.003)
[![pubmed](https://img.shields.io/badge/pubmed-24841242-blue.svg?style=flat)](https://pubmed.ncbi.nlm.nih.gov/24841242/)
[![pmc](https://img.shields.io/badge/pmc-PMC4045398-blue.svg?style=flat)](https://pmc.ncbi.nlm.nih.gov/articles/PMC4045398/)

Fraternali, F. and van Gunsteren, W.F.<br>
*An efficient mean solvation force model for use in molecular dynamics simulations of proteins in aqueous solution.*<br>
**Journal of Molecular Biology** 256 (1996) 939-948.<br>
[![doi](https://img.shields.io/badge/doi-10.1006%2Fjmbi.1996.0139-blue.svg?style=flat)](https://doi.org/10.1006/jmbi.1996.0139)
[![pubmed](https://img.shields.io/badge/pubmed-8601844-blue.svg?style=flat)](https://pubmed.ncbi.nlm.nih.gov/8601844/)

#### POPS method
Fraternali, F. and Cavallo, L.<br>
*Parameter optimized surfaces (POPS): analysis of key interactions and conformational changes in the ribosome.*<br>
**Nucleic Acids Research** 30 (2002) 2950-2960.<br>
[![doi](https://img.shields.io/badge/doi-10.1093%2Fnar%2Fgkf373-blue.svg?style=flat)](https://doi.org/10.1093/nar/gkf373)
[![pubmed](https://img.shields.io/badge/pubmed-12087181-blue.svg?style=flat)](https://pubmed.ncbi.nlm.nih.gov/12087181/)
[![pmc](https://img.shields.io/badge/pmc-PMC117037-blue.svg?style=flat)](https://pmc.ncbi.nlm.nih.gov/articles/PMC117037/)

#### POPS server
Cavallo, L., Kleinjung, J. and Fraternali, F.<br>
*POPS: A fast algorithm for solvent accessible surface areas at atomic and residue level.*<br>
**Nucleic Acids Research** 31 (2003) 3364-3366.<br>
[![doi](https://img.shields.io/badge/doi-10.1093%2Fnar%2Fgkg601-blue.svg?style=flat)](https://doi.org/10.1093/nar/gkg601)
[![pubmed](https://img.shields.io/badge/pubmed-12824328-blue.svg?style=flat)](https://pubmed.ncbi.nlm.nih.gov/12824328/)
[![pmc](https://img.shields.io/badge/pmc-PMC169007-blue.svg?style=flat)](https://pmc.ncbi.nlm.nih.gov/articles/PMC169007/)

#### POPSCOMP method and server
Kleinjung, J. and Fraternali, F.<br>
*POPSCOMP: an automated interaction analysis of biomolecular complexes.*<br>
**Nucleic Acids Research** 33 (2005) W342-W346.<br>
[![doi](https://img.shields.io/badge/doi-10.1093%2Fnar%2Fgki369-blue.svg?style=flat)](https://doi.org/10.1093/nar/gki369)
[![pubmed](https://img.shields.io/badge/pubmed-15980485-blue.svg?style=flat)](https://pubmed.ncbi.nlm.nih.gov/15980485/)
[![pmc](https://img.shields.io/badge/pmc-PMC1160130-blue.svg?style=flat)](https://pmc.ncbi.nlm.nih.gov/articles/PMC1160130/)


## Issues
Issues can be recorded on the [Issues](https://github.com/Fraternalilab/POPScomp/issues) page.


## License and Copyright
Usage of the software and server is free under the GNU General Public License v3.0.

#### Copyright Holders, Authors and Maintainers 
- 2002-2026 Franca Fraternali (author, maintainer) f.fraternali@ucl.ac.uk
- 2008-2026 Jens Kleinjung (author, maintainer) jens@jkleinj.eu

#### Contributors
- 2002 Kuang Lin and Valerie Hindie (translation to C)
- 2002 Luigi Cavallo (parametrisation)

## P.S.

#### Legacy POPScomp servers referenced in our papers
The code, servers and URL addresses in our papers have been comprehensively replaced by this GitHub repository.
For clarity, these are the legacy servers:
* mathbio.crick.mrc.ac.uk/POPSCOMP : Only for a short while the POPS and POPSCOMP servers were supported by the Francis Crick Institute.
* mathbio.nimr.mrc.ac.uk/POPSCOMP : This was the original address of the POPS and POPSCOMP servers in the Division of Mathematical Biology at the NIMR institute, which has been closed in 2015 and became a part of the Francis Crick Institute.
