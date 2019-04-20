
# POPS: Parameter OPtimised Surface of proteins and nucleic acids 
Calculation of solvent accessible surface areas (SASAs) of biopolymers,
currently proteins and nucleic acids.


## Servers
* [POPS](http://mathbio.crick.ac.uk/wiki/POPS)
* [POPSCOMP](http://mathbio.crick.ac.uk/wiki/POPSCOMP)


## References
Users publishing results obtained with the program and its applications
should acknowledge its use by the following citation:

### Implicit solvent
   Fraternali, F. and van Gunsteren, W.F.
   *An efficient mean solvation force model for use in molecular dynamics simulations of proteins in aqueous solution.*
   **Journal of Molecular Biology** 256 (1996) 939-948.

### POPS method
   Fraternali, F. and Cavallo, L.
   *Parameter optimized surfaces (POPS): analysis of key interactions and conformational changes in the ribosome.*
   **Nucleic Acids Research** 30 (2002) 2950-2960.

### POPS server
   Cavallo, L., Kleinjung, J. and Fraternali, F.
   *POPS: A fast algorithm for solvent accessible surface areas at atomic and residue level.*
   **Nucleic Acids Research** 31 (2003) 3364-3366.

### POPSCOMP server
   Kleinjung, J. and Fraternali, F.
   *POPSCOMP: an automated interaction analysis of biomolecular complexes.*
   **Nucleic Acids Research** 33 (2005) W342-W346.


## Install / Uninstall
* Dependencies
1. *libxml2* for XML parsing
2. *zlib* for compressed input files
3. *cJSON* for JSON output

* Run the 'bootstrap' shell script and follow the general 'INSTALL' instructions.

* In case the Autoconf version is too old, run
```
autoconf --version
```
and use the result to replace the version number in configure.ac
```
AC_PREREQ([2.59])
```
to reconfigure with the bootstrap script.

The configuration option '--enable-debug' creates a debuggable binary,
the option '--enable-profiling' enables profiling with 'gprof'.

For automatic code documentation execute 'doxygen doxygen.cfg'
in the 'src' directory.
Documentation files are created in 'doc/html' and 'doc/latex'.
The latex documentation is completed by executing 'make pdf' in the
'doc/latex' directory, which creates 'refman.ps' and 'refman.pdf'.

* The examples in 'tests' run by the 'make check' invocation are as follows:
- 5lff (test1 series): A peptide structure with 7 residues in one chain. 
- 1f3r (test2 series): An antibody-peptide complex with one HETATM residue.
- 1aki (test3 series): The GROMOS trajectory of a small molecule.  

## Usage
### Command Line Parameters
```
pops [--pdb ... | --pdbml ...] [OPTIONS ...]
	 INPUT OPTIONS
       Input is either a *.pdb[.gz] file (--pdb) or a *.xml[.gz] file (--pdbml).
	   --pdb <PDB input>		(mode: optional, type: char  , default: void)
	   --pdbml <PDBML input>	(mode: optional, type: char  , default: void)
	   --traj <trajectory input>	(mode: optional , type: char  , default: void)
	   --zipped                 (mode: opional , type: bool  , default: false)
	 MODE OPTIONS
	   --coarse			(mode: optional , type: no_arg, default: off)
	   --hydrogens			(mode: optional , type: no_arg, default: off)
	   --multiModel			(mode: optional , type: no_arg, default: off)
	   --partOcc			(mode: optional , type: no_arg, default: off)
	   --rProbe <probe radius [A]>	(mode: optional , type: float , default: 1.4)
	   --silent			(mode: optional , type: no_arg, default: off)
	 OUTPUT OPTIONS
	   --popsOut <POPS output>	(mode: optional , type: char  , default: pops.out)
	   --popstrajOut <POPS output>	(mode: optional , type: char  , default: popstraj.out)
	   --popsbOut <POPSb output>	(mode: optional , type: char  , default: popsb.out)
	   --popsbtrajOut <POPSb output>(mode: optional , type: char  , default: popsbtraj.out)
	   --sigmaOut <SFE output>	(mode: optional , type: char  , default: sigma.out)
	   --sigmatrajOut <SFE output>	(mode: optional , type: char  , default: sigmatraj.out)
	   --interfaceOut		(mode: optional , type: no_arg, default: off)
	   --compositionOut		(mode: optional , type: no_arg, default: off)
	   --typeOut			(mode: optional , type: no_arg, default: off)
	   --topologyOut		(mode: optional , type: no_arg, default: off)
	   --atomOut			(mode: optional , type: no_arg, default: off)
	   --residueOut			(mode: optional , type: no_arg, default: off)
	   --chainOut			(mode: optional , type: no_arg, default: off)
	   --neighbourOut		(mode: optional , type: no_arg, default: off)
	   --parameterOut		(mode: optional , type: no_arg, default: off)
	   --noTotalOut			(mode: optional , type: no_arg, default: off)
	   --noHeaderOut		(mode: optional , type: no_arg, default: off)
	   --padding			(mode: optional , type: no_arg, default: off)
	   --rout			(mode: optional , type: no_arg, default: off)
	   --jsonOut			(mode: optional , type: no_arg, default: off)
	 INFO OPTIONS
	   --cite			(mode: optional , type: no_arg, default: off)
	   --version			(mode: optional , type: no_arg, default: off)
	   --help
```

### Short Description of Command Line Parameters
* pdb : input format is the classical PDB format
* pdbml : input format is the XML format of the PDB database
* trajInFileName : trajectory input file
* zipped : the input file is compressed with gzip or similar (zlib compatible)
* coarse : Calpha-only computation [0,1]
* hydrogens : hydrogens [0,1]
* multiModel : input with multiple models
* rProbe : probe radius (in Angstrom) */
* silent : suppress stdout
* sasaOutFileName : output of SASA values for single structure
* sasatrajOutFileName : output of SASA values for trajectory
* bsasaOutFileName : output of buried SASA for single structure
* bsasatrajOutFileName : output of buried SASA for 
* compositionOut : print molecule composition
* sigmaOutFileName : output of sigma values for single structure 
* sigmatrajOutFileName : output of sigma value for trajectory
* interfaceOutFileName : pairlist of contact residues between distinct chains, where only the contact residue with the shortest distance is listed
* interfaceOut : output of interface residue pairs
* typeOut : print atom/residue types
* topologyOut : print molecule topology
* atomOut : print atom areas
* residueOut : print residue areas
* chainOut : print chain areas
* neighbourOutFileName : output of neighbour list 
* neighbourOut : print neighbour list 
* parameterOutFileName : output of POPS parameters
* parameterOut : print atom parameters
* noTotalOut : suppress output of total area (for benchmarking)
* noHeaderOut : suppress output headers (for benchmarking)
* padding : add lines to pad the missing hydrogen atom lines (for benchmarking)
* rout : output of R-compatible output 


## Input file formats
POPS accepts two formats of the PDB database:
the classical PDB format (*.pdb) using input parameter '--pdb' and the
newer XML format (*.xml) using input parameter '--pdbml';
the latter is the preferred option for forward compatibility.
Structures are usually stored in compressed (zipped) form and
POPS can read the compressed files directly.
POPS depends therefore on the libxml2 and lzip libraries,
whose presence is checked at the configuration stage.


## Program Design
* PDB input file (read input) 
* Groups (assign atom group id)
* Types (assign GROMOS atom type id)
* Topology (compute molecular topology)
* SASA (compute SASA values at atoms, chain and molecule level)
* Output (print output)
* From version 3.0 the POPSCOMP method has been built into the POPS program
    so that complex interfaces are computed automatically.


### Determine the atom and residue types according to 'sasa_data.h' ('type'),
Possible complications in the handling of PDB files and the way it
is dealt with by POPSc are listed here.
* Hydrogen atoms : Hydrogen atoms are skipped.
* Alternative locations : Alternative locations are skipped.
* Multiple chains : The chain number of each atom is printed in the 'atomOut' output.
  Using the option 'chainOut' prints the chain areas.
* Multiple models : Only the first model is used.
* Missing atoms/residues : No detection of missing atoms/residues.
* Unrecognised standard atoms/residues (ATOM entries): 
  The program issues a warning and/or exits.
  Some programs create non-standard atom entries, which are correct within
  the limits of the PDB format definition, but differ from the standard
  format that is used by the PDB for amino acids and nucleotides.
  POPS expects the standard format as in these examples:
* amino acid
    ATOM     41  CA  PRO A  69      -3.665   2.182  -6.278  1.00  0.00           C
* nucleotide
    ATOM     33  P     G A  17       4.130  21.140  -1.230  1.00  0.00           P 
* Hetero atoms (HETATM entries): These are treated differently from the
  standard ATOM entries. Atom names are matched against the standard
  atom names [' N  ','  CA',' C  ',' O  '] to recognise modified
  amino acids, otherwise atom names are matched against the elements
  ['C','N','O','P','S'], in which case these atoms are renamed to
  ['C_','N_','O_','P_','S_'] to indicate their generic nature.
  Unrecognised HETATM entries are skipped!
* For a complete list of residues see:
    [monomers](ftp://ftp.wwpdb.org/pub/pdb/data/monomers/)

### Atomgroup assignment (code file 'atomgroup').
POPS* uses a specific sigma parameter for each atom type.
Atom types with the same sigma can be grouped together, 
yielding a group-specific sigma parameter.

### Type assignment for atoms and residues (code file 'type').
The data structure that defines the atom and residue types is specified
as data type 'ConstantSasa' in 'sasa_const.h', 
while the actual parametrisation is kept in 'sasa_data.h'.
'sasa_data.h' defines an array 'ConstantSasa constant_sasa_data[]'.
Each element of this array is a complete POPS parametrisation.
At the moment there are two parametrisations, element '0' is the
atomic level, element '1' is the coarse-grained residue level.
The pointer 'constant_sasa' in the main program ('pops') points
to either element, which defines the parametrisation to be used.
Other parametrisations can be defined by appending new elements
to the 'constant_sasa_data[]' array.

### Determine the molecular Topol (code file 'topol').
The determination of the molecular topology proceeds progressively
from 1-2 interactions (bonds) over 1-3 interactions (angles) to
1-4 interactions (torsions). Any other pairwise atom pair relationship 
is by definition non-bonded. The 1-(i+1) interactions are assembled
from all pariwise combinations of 1-i interaction by testing identical
atom overlaps. The determination of ring topologies is a 
painfully slow and error-prone procedure if not implemented ingeniously.
Therefore, the ring status of each atom (0=no ring, 1=ring) has been
included in the 'ConstantSasa' data type so that there is no need to
determine it. This means the treatment of cyclic HETATM molecules
requires specific attention.

### SASA: Compute AtomSasa, ResSasa, ChainSasa and MolSasa (code file 'sasa').
The SASA calculation has been broken down into elementary equations (as
inline code) and a progressive strategy to approximate atom SASAs.
Starting from the SASA of the isolated atom, overlap areas are subtracted
progressively for all 1-(i+n) interactions. Atomic SASAs are summed up
to yield residue, chain and molecule SASAs.

### Output
Flags controlling the amount and type of output.
The default ouput is the overall SASA of the molecule, including the
hydrophilic and hydrophobic contributions. More output can be produced
with the '--*Out' command line flags.

### ATOM SASAs
* AtomNr : atom number in molecular coordinate file
* AtomNe : atom name in molecular coordinate file
* ResidNe : residue name in molecular coordinate file
* Chain : chain name in molecular coordinate file (can be void)
* ResidNr : residue number in molecular coordinate file
* SASA/A^2 : solvent accessible surace area in Angstrom^2 units
* Q(SASA) : quotient of SASA and Surf (below), i.e. the fraction of SASA
* N(overl) : number of overlaps with atom neighbours
* AtomTp : atom type code (GROMOS van der Waals atom type)
* AtomGp : atom group code: positive=1, negative=2, polar=3, aromatic=4, aliphatic=5
* Surf/A^2 : surface area of isolated atom

### RESIDUE SASAs
* ResidNe : residue name in molecular coordinate file
* Chain : chain name in molecular coordinate file (can be void)
* ResidNr : residue number in molecular coordinate file
* Phob/A^2 : hydrophobic solvent accessible surace area in Angstrom^2 units
* Phil/A^2 : hydrophilic solvent accessible surace area in Angstrom^2 units
* Total/A^2 : total solvent accessible surace area in Angstrom^2 units
* Q(SASA) : quotient of SASA and Surf (below), i.e. the fraction of SASA
* N(overl) : number of overlaps with residue neighbours
* Surf/A^2 : surface area of isolated residue

### CHAIN SASAs
* Chain : chain number 
* Id : chain name in molecular coordinate file (can be void)
* AtomRange : range of atom numbers in chain
* ResidRange : range of residue numbers in chain
* Phob/A^2 : hydrophobic solvent accessible surace area in Angstrom^2 units
* Phil/A^2 : hydrophilic solvent accessible surace area in Angstrom^2 units
* Total/A^2 : total solvent accessible surace area in Angstrom^2 units

### MOLECULE SASAs
* Phob/A^2 : hydrophobic solvent accessible surace area in Angstrom^2 units
* Phil/A^2 : hydrophilic solvent accessible surace area in Angstrom^2 units
* Total/A^2 : total solvent accessible surace area in Angstrom^2 units

### bSASAs
The same type of information as for the SASAs is produced for bSASAs,
which are the buried SASAs of side-chain atoms only.
bSASAs serve as a component in the definition of a structural 
environment (3D profile).
* Phob/A^2 : solvent accessible surface area buried by hydrophobic neighbour atoms in Angstrom^2 units
* Phil/A^2 : solvent accessible surface area buried by hydrophilic neighbour atoms in Angstrom^2 units
* Total/A^2 : solvent accessible surface area buried by all neighbour atoms in Angstrom^2 units

### INTERFACE
Information of interface residues in complexes can be obtained with switch
--interfaceOut. For each atom of the input structure, the nearest neighbour
is computed based on the neighbour list and printed if the two contact
atoms are in different chains. Note that this pairlist is not symmetric:
atom distance A-B might be shortest for A, but distance B-C might be shortest
for B. Also, the list might change with protein motion. Therefore it depends
on the interpretation to be derived from the interface whether it is suitable
to use the atom pairs, the residue pairs and/or their distances.

### Background Information
The area equation is defined by a product ∏ of terms that estimate
the reduction of SASA of atom i by the overlap with its 
neighbours j (Hasel at al., 1988):
```
∏^N_i=1 [1 - (p_i p_ij b_ij (r_ij) / S_i)].
```
* i is the atom for which the POPS area is computed, j is any of N neighbour atoms.
* p_i  is an atom type specific SASA parameter.
* p_ij is a sphere overlap parameter depending on the degree of bonding between i and j (1-2, 1-3, 1-4, 1-5, ...).
* b_ij is a geometric construct based on the radii and distance (r_ij) of i and j.
* S_i is the SASA of the free atom i (no neighbours).

The atom specific parameters (radii, SASAs) are listed on the web server's
'Parameters' files for the atoms of all standard PDB residues, followed by the connectivity parameters (p_ij, b_ij)
and the solvent radius for water.


### Return values
* 0 : Clean termination.
* 1 : Error:
  - Input error (wrong file name).
  - Too few arguments.
  - Unknown standard residue or atom. 
  - Too short atom distance.
  - Less than 2 atoms, 2 bonds, 2 angles or 2 torsions in input structure.
* Assertions: Internal data flow is partially controlled by assertions.
  The 'assert' macro does not return an error value itself,
  but the program exits with a (non-zero) error code and a comment.
* Warning messages: Error and warninig output is directed to 'stderr'.
* Progress messsages: Progress output is sparse and directed to 'stdout'.
* Default result output: By default the program prints results to 'pops.out'.


## Code
* [POPS code](https://github.com/Fraternalilab/POPS)
* [Latest Release](https://github.com/Fraternalilab/POPS/releases/latest)

### Github branches
* master : the distributed version of POPS
* PDBML : parses XML format instead of PDB format
* json : prints POPS output in JSON format
* popscomp : computes results of POPSCOMP protocol
* Rout : prints POPS output in R compatible table format


## Contact
* franca.fraternali@kcl.ac.uk
* jens@jkleinj.eu


## Copyright
* 2002-2018 Franca Fraternali (program author)
* 2008-2018 Jens Kleinjung (modular C code)
* 2002 Kuang Lin and Valerie Hindie (translation to C)
* 2002 Luigi Cavallo (parametrisation)


## Availability
The program is made available under the GNU Public License for academic
scientific purposes, under the condition that proper acknowledgement
is made to the authors of the program in publications resulting from the use
of the program.


## License
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

