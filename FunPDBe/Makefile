#===============================================================================
# Makefile : POPScomp pipeline for FunPDBe 
# Copyright (C) 2019 Jens Kleinjung
#===============================================================================

#_______________________________________________________________________________
## start-up macros
## root directory
ROOTDIR = /home/jkleinj
## POPSC program
POPSC = ${ROOTDIR}/0/POPSCOMP/POPSC/src/pops
## database directory
DBDIR = ${ROOTDIR}/database
## home of PDB database
PDBDIR = ${DBDIR}/XML
## PDB subdirectories
PDBSUBDIR = $(notdir $(shell find ${PDBDIR}/ -type d))
## home of JSON output files
JSONDIR = ${DBDIR}/JSONTEST
## JSON subdirectories
JSONSUBDIR = $(notdir $(shell find ${JSONDIR}/ -type d))

#_______________________________________________________________________________
## suffix rules
.SUFFIXES :

.SUFFIXES : .gz .json

%.json : %.gz
	$(info Running POPS over ${PDBIDS})
	${POPSC} --pdbml $< --outDirName --zipped --jsonOut ${JSONDIR}/$X --silent > /dev/null 2> popsPDBML.e

#_______________________________________________________________________________
## targets
all : pdbupdate jsondir jsonsubdir mkpops

pdbupdate :
	rsync -rlpt -v -z --delete rsync.ebi.ac.uk::pub/databases/pdb/data/structures/divided/XML ${DBDIR}/

jsondir :
	$(info --- Creating output directory.)
	mkdir -p ${JSONDIR}

jsonsubdir :
	$(info --- Creating output sub-directories.)
	cd ${JSONDIR} ; \
	for X in ${PDBSUBDIR} ; do \
		mkdir -p ${JSONDIR}/$$X ; done

mkpops : 
	cd ${PDBDIR}; \
	for X in ${PDBSUBDIR} ; do \
		${PDBIDS} = $(patsubst %.xml.gz,%.json,$(wildcard *.gz)) ; \
		$(info ${PDBIDS}) ; \
		${MAKE} run_pops "PDBID = ${PDBID}" \
						 "X = $$X" \
						 "PDBIDS = ${PDBIDS}" ; \
    done

run_pops : ${PDBID}/$X/${PDBIDS}


#===============================================================================

