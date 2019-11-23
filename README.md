## Readme
The POPS program computes the Solvent Accessible Surface Area (SASA)
of a given PDB structure. If the structure is composed of more than one chain
containing protein or RNA/DNA, POPScomp creates internally all pair combinations
of chains to compute the buried SASA upon complexation. Details of those routines
are explained in the published papers on POPS and POPSCOMP.

POPScomp (*Shiny* app) shows tabs for atom, residue, chain and molecule SASAs.
The tables are initialised without any values and therefore the user sees
the table header and below the notice 'No data available in table'.
After selecting a PDB identifier or file and pressing 'run POPScomp',
the sever runs the POPS program on components of the PDB file
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

Results will be stored on the server for maximally one day.
For permanent storage, please download your results via the
'Download' buttons.


## About
This is version 3.0.0 of the [POPScomp server](http://popscom.org:3838).
The server automatically recognises PDB identifiers and multi-chain structures.
Output comprises downloadable SASA tables and graphs shown on the Shiny server pages.

### Packages
Since POPScomp 3.0 (04.2019), the packages *POPSC*, *POPSR*, *FunPDBe* and
*Docker* are part of the POPScomp project.
1. *POPSC*: A GNU Autotools package of the POPS C program that computes SASA
  for a given structure.
2. *POPSR*: An R package containing an R program that
    - splits complexes into single and pair components to compute buried SASA
	  using POPSC and 
    - provides a Shiny server as interface to the POPScomp (C and R) programs.
3. *FunPDBe*
Scripts to run POPScomp over the PDB database and feed the output into
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
docker push eu.gcr.io/high-hue-217311/popscomp:0.1
```

### POPScomp server
* [POPScomp](http://popscomp.org:3838)

### FunPDBe
POPScomp is part of the [FunPDBe resources](https://www.ebi.ac.uk/pdbe/funpdbe/deposition).

### License
Usage of the software and server is free, the code license is GPL3.

### Copyright Holders, Authors and Maintainers 
- 2002-2019 Franca Fraternali (franca.fraternali@kcl.ac.uk)
- 2008-2019 Jens Kleinjung (jens@jkleinj.eu)

### Contributors
- 2002 Kuang Lin and Valerie Hindie (translation to C)
- 2002 Luigi Cavallo (parametrisation)

