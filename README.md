## Readme
The POPS program computes the Solvent Accessible Surface Area (SASA)
of a given PDB structure. If the structure is composed of more than one chain
containing protein or RNA/DNA, POPScomp creates internally all pair combinations
of chains to compute the buried SASA upon complexation. Details of those routines
are explained in the published papers on POPS and POPSCOMP.

POPScomp shows tabs for atom, residue, chain and molecule SASAs.
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
Upon high demand the storage time might be reduced 30 minutes.
Please download your results via the 'Download' buttons.


## About
This is version 3.0.0 of the [POPScomp server](href="http://popscom.org:3838).
The server automatically recognises PDB identifiers and multi-chain structures.
Output comprises downloadable SASA tables and graphs shown on the Shiny server pages.

The POPScomp server is based on two software packages:
1. A GNU Autotools package of the POPS C program that computes SASA for a given structure.
2. An R package containing an R program that a) splits complexes into single and pair
  components to compute buried SASA using POPSC and b) provides a Shiny server to interface the R program.
Since April 2019 the POPS program (POPSC) and the POPScomp Shiny server (POPSR)
are being co-developed.

Source code and detailed information can be found on Fraternali lab's
[POPScomp GitHub page](href="https://github.com/Fraternalilab/POPScomp).
Please use that site for bug reports and other comments.

POPScomp is part of the [FunPDBe resources](href="https://www.ebi.ac.uk/pdbe/funpdbe/deposition).

Usage of the server is free, the code license is GPL3.

Authors:
- Franca Fraternali (franca.fraternali@kcl.ac.uk)
- Jens Kleinjung (jens@jkleinj.eu)
        )
