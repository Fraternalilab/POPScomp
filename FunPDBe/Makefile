#===============================================================================
# Makefile : POPScomp pipeline for FunPDBe 
# Invokes Makefile.subdir for each subdirectory of PDBML.
# There are three 'mkpops' targets to improve parallelisation, because the file
# list per subdirectory lingers obviously on the slowest one in regard to POPS
# execution and that means only one CPU is used during that time.
# Copyright (C) 2019-2020 Jens Kleinjung
#===============================================================================

#_______________________________________________________________________________
## start-up macros
## root directory
ROOTDIR = /home/jkleinj
## POPSC program
POPSC = ${ROOTDIR}/0/POPScomp/POPSC/src/pops
POPSFLAGS = --zipped --jsonOut --silent
## database directory
DBDIR = ${ROOTDIR}/0/databases
## home of PDB database
PDBDIR = ${DBDIR}/XML
## PDB subdirectories
PDBSUBDIR = $(notdir $(shell find ${PDBDIR}/ -type d))
PDBSUBDIR1 = $(notdir $(shell find ${PDBDIR}/[a-h]? -type d))
PDBSUBDIR2 = $(notdir $(shell find ${PDBDIR}/[i-p]? -type d))
PDBSUBDIR3 = $(notdir $(shell find ${PDBDIR}/[q-z]? -type d))
## home of JSON output files
JSONDIR = ${DBDIR}/JSON
## JSON subdirectories
JSONSUBDIR = $(notdir $(shell find ${JSONDIR}/ -type d))
## home of JSON validation output files
JSONVALDIR = ${DBDIR}/JSONVAL
## JSON validation subdirectories
JSONVALSUBDIR = $(notdir $(shell find ${JSONVALDIR}/ -type d))
## FunPDBe directory
FUNPDBEDIR = ${ROOTDIR}/0/POPScomp/FunPDBe

#_______________________________________________________________________________
## targets
## pdbupdate: Synchronise local PDB database with EBI's.
## jsondir: Create JSON output directory if not existing.
## jsonsubdir: Create JSON (output) subdirectory structure, mirroring the PDB's.
## linkinput: Create symbolic links to all (input) PDBML entries,
##              required to derive the dependency/target file names for Make.
## mkpops: Run POPS over all subdirectories and PDBML entries, only if
##           JSON output is missing or older than PDBML input.
## rename: Convert JSON output filenames to lower case, otherwise Make would not
##           know from the suffix rule that the JSON outtput file exists.
## validate: Validate the JSON output against the FunPDBe schema
##             (funpdbe_client.py).
##
## Note: PDBSUBDIR is also used for JSONDIR, because it is the reference array
##         that is completely known from the start.

all : pdbupdate  jsondir jsonsubdir linkinput mkpops rename jsonvaldir jsonvalsubdir linkoutput validate rsyncval mkupload
.PHONY : pdbupdate jsondir jsonsubdir linkinput mkpops mkpops1 mkpops2 mkpops3 rename jsonvaldir jsonvalsubdir linkoutput validate rsyncval mkupload


pdbupdate :
	rsync -rlpt -v -z --delete rsync.ebi.ac.uk::pub/databases/pdb/data/structures/divided/XML/ ${DBDIR}/


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
	cd ${PDBDIR} ; find ?? -depth -name *.xml.gz -mtime -28 | xargs -i bash -c "ln -s ${PDBDIR}/{} ${JSONDIR}/{}"


#.PHONY: mkpops $(PDBSUBDIRS)

#mkpops : $(PDBSUBDIRS)

#$(PDBSUBDIRS) :
#	${MAKE} -f ../Makefile.subdir all "ROOTDIR = ${ROOTDIR}" "POPSC = ${POPSC}" "POPSFLAGS = ${POPSFLAGS}" "PDBDIR = ${PDBDIR}" $@ 

mkpops :
	$(info --- Run POPS on all structures.)
	for X in ${PDBSUBDIR} ; do \
		cd ${JSONDIR}/$$X ; \
		${MAKE} -i -f ../Makefile.subdir all "ROOTDIR=${ROOTDIR}" "POPSC=${POPSC}" "POPSFLAGS=${POPSFLAGS}" "PDBDIR=${PDBDIR}" "X=$$X" ; \
	done

mkpops1 :
	$(info --- Run POPS on structures [a-h]?.)
	for X in ${PDBSUBDIR1} ; do \
		cd ${JSONDIR}/$$X ; \
		${MAKE} -i -f ../Makefile.subdir all "ROOTDIR=${ROOTDIR}" "POPSC=${POPSC}" "POPSFLAGS=${POPSFLAGS}" "PDBDIR=${PDBDIR}" "X=$$X" ; \
	done

mkpops2 : 
	$(info --- Run POPS on structures [i-p]?.)
	for X in ${PDBSUBDIR2} ; do \
		cd ${JSONDIR}/$$X ; \
		${MAKE} -i -f ../Makefile.subdir all "ROOTDIR=${ROOTDIR}" "POPSC=${POPSC}" "POPSFLAGS=${POPSFLAGS}" "PDBDIR=${PDBDIR}" "X=$$X" ; \
	done

mkpops3 : 
	$(info --- Run POPS on all structures [q-z]?.)
	for X in ${PDBSUBDIR3} ; do \
		cd ${JSONDIR}/$$X ; \
		${MAKE} -i -f ../Makefile.subdir all "ROOTDIR=${ROOTDIR}" "POPSC=${POPSC}" "POPSFLAGS=${POPSFLAGS}" "PDBDIR=${PDBDIR}" "X=$$X" ; \
	done

rename :
	cd ${JSONDIR} ; find ?? -depth -name "*[A-Z]*".json -mtime -28 | xargs -n 1 rename -f -v 'y/A-Z/a-z/' {} 

#_______________________________________________________________________________
## validation
jsonvaldir :
	$(info --- Creating validation output directory.)
	mkdir -p ${JSONVALDIR}

jsonvalsubdir :
	$(info --- Creating validation output sub-directories.)
	cd ${JSONVALDIR} ; \
	for X in ${PDBSUBDIR} ; do \
		mkdir -p ${JSONVALDIR}/$$X ; \
	done

linkoutput :
	$(info --- Linking POPS output files to validation directory.)
	cd ${JSONDIR} ; find ?? -depth -name *.json -mtime -28 | xargs -i bash -c "ln -s ${JSONDIR}/{} ${JSONVALDIR}/{}"


validate :
	$(info --- Validating all output files.)
	for X in ${PDBSUBDIR} ; do \
		cd ${JSONVALDIR}/$$X ; \
		${MAKE} -i -f ../Makefile.validate validate "JSONDIR=${JSONDIR}" "X=$$X" "FUNPDBEDIR=${FUNPDBEDIR}" ; \
	done

rsyncval :
	$(info --- Get validated targets via rsync from JSONVAL on mito.)
	for X in ${PDBSUBDIR} ; do \
		cd ${JSONVALDIR}/$$X ; \
		${MAKE} -i -f ../Makefile.validate rsyncval "JSONDIR=${JSONDIR}" "X=$$X" "FUNPDBEDIR=${FUNPDBEDIR}" ; \
	done

#_______________________________________________________________________________
mkupload :
	$(info --- Uploading validated output files.)
	${MAKE} -i -f Makefile.upload upload "JSONVALDIR = ${JSONVALDIR}"


#===============================================================================

