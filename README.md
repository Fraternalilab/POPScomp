# POPScomp
POPScomp is a program for the computation of Solvent Accessible Surface Area (SASA)
of proteins, DNA/RNA and their complexes.
The program is composed of two components, POPSC and POPSR, which are being
developed together:
- POPSC is a C program that computes SASA for a given structure.
- POPSR is an R program that splits complexes into single and pair components
  to compute buried SASA using POPSC.

## Server
We provide a server at [popscomp.org](http://popscomp.org:3838)
that automatically recognises PDB identifiers and multi-chain structures,
creates various output tables and graphs.

## Code
The source code is available as GitHub repository of the Fraternali lab:
[POPScomp repository](https://github.com/Fraternalilab/POPScomp).
Please use that site for bug reports and other comments.

## License
Usage of the server is free, the code license is GPL3.

## Copyright
(C) 2019 Jens Kleinjung and Franca Fraternali

