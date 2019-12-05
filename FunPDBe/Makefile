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
JSONDIR = ${DBDIR}/JSON
## JSON subdirectories
JSONSUBDIR = $(notdir $(shell find ${JSONDIR}/ -type d))
## selected subdirectory, will become meaningful during computation
SELSUBDIR = 
## PDB IDs in selected subdirectory, i.e. those to process now
PDBIDS = $(notdir $(shell find ${PDBDIR}/${SELSUBDIR} -type f))


#_______________________________________________________________________________
## suffix rules
.SUFFIXES :

.SUFFIXES : .gz .json

##%/${PDBIDS}.json : %/${PDBIDS}.xml.gz
##    $(info Running POPS over ${PDBIDS})
##	cd ${DATABASE}; ${POPSC} --pdbml $< --outDirName --zipped --jsonOut --silent > /dev/null 2> popsPDBML.e

## PDB files
#PDB = $(shell find ${database}/XML/ -type f -name *.xml.gz)

#_______________________________________________________________________________
## targets
all : pdbupdate jsondir jsonsubdir runpops

pdbupdate :
	cd ${DATABASE}; \
	rsync -rlpt -v -z --delete rsync.ebi.ac.uk::pub/databases/pdb/data/structures/divided/XML ${DBDIR}/

jsondir :
	$(info --- Creating output directory.)
	mkdir -p ${JSONDIR}

jsonsubdir :
	$(info Creating output sub-directories.)
	for X in ${PDBSUBDIR} ; do \
		mkdir -p ${JSONDIR}/$$X ; done

runpops : ${PDBIDS}

${PDBIDS} : 
    $(info --- Running POPS.)
	##$(info MUTS = ${MUTS})
    ##cd ${EXP}; \
	##${POPSC} --pdbml $@ --outDirName --zipped --jsonOut --silent > /dev/null 2> popsPDBML.e

#===============================================================================

