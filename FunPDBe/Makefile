#===============================================================================
# Makefile : POPScomp pipeline for FunPDBe 
# Invokes Makefile.subdir for each subdirectory
# Copyright (C) 2019 Jens Kleinjung
#===============================================================================

#_______________________________________________________________________________
## start-up macros
## root directory
ROOTDIR = /home/jkleinj
## POPSC program
POPSC = ${ROOTDIR}/0/POPScomp/POPSC/src/pops
POPSFLAGS = --zipped --jsonOut --silent
## database directory
DBDIR = ${ROOTDIR}/database
## home of PDB database
PDBDIR = ${DBDIR}/XML
## PDB subdirectories
#PDBSUBDIR = $(notdir $(shell find ${PDBDIR}/ -type d))
PDBSUBDIR1 = $(notdir $(shell find ${PDBDIR}/[a-h]? -type d))
PDBSUBDIR2 = $(notdir $(shell find ${PDBDIR}/[i-p]? -type d))
PDBSUBDIR3 = $(notdir $(shell find ${PDBDIR}/[q-z]? -type d))
## home of JSON output files
JSONDIR = ${DBDIR}/JSON
## JSON subdirectories
JSONSUBDIR = $(notdir $(shell find ${JSONDIR}/ -type d))

#_______________________________________________________________________________
## targets
## pdbupdate: Synchronise local PDB database with EBI
## jsondir: Create JSON output directory if not existing
## jsonsubdir: Create JSON (output) subdirectory structure, mirroring the PDB structure
## linkinput: Create symbolic links to all (input) PDBML entries
## mkpops: Run POPS over all subdirectories and PDBML entries, only if JSON output
##           is null or older than PDBML input.
## rename: Convert JSON output filenames to lower case, otherwise Make would not
##           know from the suffix rule that the JSON outtput file exists
## Note: PDSUBDIR is also used for JSONDIR, because it is the reference array
##         that is completely known from the start

all : pdbupdate jsondir jsonsubdir getstr mkpops1 mkpops2 mkpops3 rename

pdbupdate :
	rsync -rlpt -v -z --delete rsync.ebi.ac.uk::pub/databases/pdb/data/structures/divided/XML ${DBDIR}/

jsondir :
	$(info --- Creating output directory.)
	mkdir -p ${JSONDIR}

jsonsubdir :
	$(info --- Creating output sub-directories.)
	cd ${JSONDIR} ; \
	for X in ${PDBSUBDIR} ; do \
		mkdir -p ${JSONDIR}/$$X ; \
	done

linkinput :
	$(info --- Linking structure files to output directory.)
	cd ${PDBDIR} ; find ?? -depth -name *.xml.gz | xargs -i bash -c "ln -s ${PDBDIR}/{} ${JSONDIR}/{}"

mkpops1 : 
	$(info --- Run POPS on all structures.)
	for X in ${PDBSUBDIR1} ; do \
		cd ${JSONDIR}/$$X ; \
		${MAKE} -f ../Makefile.subdir all "ROOTDIR = ${ROOTDIR}" "POPSC = ${POPSC}" "POPSFLAGS = ${POPSFLAGS}"; \
	done

mkpops2 : 
	$(info --- Run POPS on all structures.)
	for X in ${PDBSUBDIR2} ; do \
		cd ${JSONDIR}/$$X ; \
		${MAKE} -f ../Makefile.subdir all "ROOTDIR = ${ROOTDIR}" "POPSC = ${POPSC}" "POPSFLAGS = ${POPSFLAGS}"; \
	done

mkpops3 : 
	$(info --- Run POPS on all structures.)
	for X in ${PDBSUBDIR3} ; do \
		cd ${JSONDIR}/$$X ; \
		${MAKE} -f ../Makefile.subdir all "ROOTDIR = ${ROOTDIR}" "POPSC = ${POPSC}" "POPSFLAGS = ${POPSFLAGS}"; \
	done

rename :
	cd ${JSONDIR} ; find ?? -depth -name *[A-Z]*.json | xargs -n 1 rename -v 'y/A-Z/a-z/' {}

#===============================================================================


