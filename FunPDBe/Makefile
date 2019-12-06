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

#_______________________________________________________________________________
## suffix rules
.SUFFIXES :

.SUFFIXES : .gz .json

%.json : %.gz
	$(info Running POPS over ${PDBIDS})
	${POPSC} --pdbml $< --outDirName --zipped --jsonOut ${JSONDIR}/$X --silent > /dev/null 2> popsPDBML.e

#_______________________________________________________________________________
## targets
## pdbupdate: Synchronise local PDB database with EBI
## jsondir: Create JSON output directory if not existing
## jsonsubdir: Create JSON subdirectory structure, mirroring the PDB structure
## mkpops: Run POPS over all subdirectories and PDB entries without a JSON output file
## rename: Convert JSON output filenames to lower case, otherwise Make would not
##           know from the suffix rule that the JSON outtput file exists
##         Also add a symbolic link to the source structure on the PDB database
all : pdbupdate jsondir jsonsubdir mkpops rename

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

mkpops : 
	$(info --- Run POPS on all structures.)
	cd ${JSONDIR} ; \
	for X in ${PDBSUBDIR} ; do \
		ls *.json ; \
		${PDBIDS} = $(patsubst %.xml.gz,%.json,$(wildcard *.gz)) ; \
	done

#	for X in ${PDBSUBDIR} ; do \
#		${PDBIDS} = $(patsubst %.xml.gz,%.json,$(wildcard *.gz)) ; \
#		${MAKE} run_pops "PDBID = ${PDBID}" "X = $$X" "PDBIDS = ${PDBIDS}" ; \
#	done

run_pops : ${PDBDIR}/$X/${PDBIDS}


rename :
##	cd ${JSONDIR} ; find ?? -depth -name *.json | xargs -n 1 rename -v 's/(.*)\/([^\/]*)/$1\/\L$2/' {} ; 
	cd ${JSONDIR} ; find ?? -depth -name *.json | cut -f 1 -d '.' | xargs -i bash -c "ln -s ${PDBDIR}/{}.xml.gz ${JSONDIR}/{}.xml.gz"


#===============================================================================

