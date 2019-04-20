# *POPScomp*
To make the text below more readable, programs (in the sense of
data processing routines) are formatted in **bold**, while interfaces
to their input/output streams are denoted in *italic*.

The *POPScomp* package provides the R Shiny interface to **POPS**,
which is a C program and neeeds to be installed and compiled separately;
the latter should be present under 'POPScomp/bin/pops' as a link.
See 'POPScomp/bin/README.md' for more information.

The **POPS** program computes the Solvent Accessible Surface Area (SASA)
of a given PDB structure. If the structure was composed of more than one chain
containing protein or RNA/DNA, *POPScomp* creates internally all pair combinations
of chains to compute the buried SASA upon complexation. Details of the procedure
are explained in the published papers on POPS and POPSCOMP.
In that regard, *POPScomp* (starting with *POPS* version 3.0) replaces both
the original *POPS* and *POPSCOMP* web server interfaces (written in PHP) as well
as the **POPSCOMP** script for multiple chain processing (written in Perl).

*POPScomp* shows several tabs for atom, residue, chain and molecule SASAs.
The tables are initialised without any values and therefore the user sees
the table header and below the notice "Showing 0 to 0 of 0 entries".
After selecting a PDB file and launching "run POPScomp", the sever runs
the **POPS** program on components of the PDB file. Because that computation
is a system call, the *POPScomp* server sees only the exit code of that system call;
the exit code will be shown below the "run POPScomp" button:

* 0 - Success
* 1 - Catchall for general errors
* 2 - Misuse of shell builtins (according to Bash documentation)
* 126 - Command invoked cannot execute
* 127 - “command not found”
* 128 - Invalid argument to exit
* 128+n - Fatal error signal “n”
* 130 - Script terminated by Control-C
* 255\* - Exit status out of range 

