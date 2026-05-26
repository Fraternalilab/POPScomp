#!/usr/bin/env Rscript

#===============================================================================
# POPSR package
# popscomp.R: Implementation of the POPSCOMP functionality,
# i.e. processing of complex structures to compute SASA difference values.
# Returns a list of POPS output files for single-chain and pair-chain structures
#   plus a list of buried SASA values.
#
# (C) 2019-2026 Jens Kleinjung and Franca Fraternali
#===============================================================================

library("bio3d")
## For atom selection mechanisms see:
## http://thegrantlab.org/bio3d/tutorials/structure-analysis
## The 'pdbsplit' function used below to split the structure into chains
##   works only on PDB files (not MMCIF).
##   For MMCIF files, chain splitting is performed here via 'gemmi'
##   on the command line.
## Note that 'pops' (POPSC) is also capable of reading MMCIF files
##   via the in-built 'gemmi' library, but chain splitting is implemented
##   only in this 'popscomp' (POPSR) part of the program suite.

library("optparse")

#_______________________________________________________________________________
## POPScomp function implemented in R
## The following prefixes are used to label the sections and output files
##   for clarity (DIFF output files are called 'delta' for historic reasons):
## ID: the default '--popsr' prefix of POPS for the unmodified input PDB
##      (computed by 'input$popscomp' function in 'app.R')
## ISO: POPS on isolated chains
## PAIR: POPS on paired chains
## DIFF: difference between sum of isolated chain SASA and paired chain SASA

option_list = list(
  make_option(c("-p", "--pdb"), type = "character", default = NULL,
              help = "local PDB file (upload)", metavar = "character"),
  make_option(c("-i", "--id"), type = "character", default = NULL,
              help = "PDB identifier (download)", metavar = "character"),
  make_option(c("-m", "--mmcif"), type = "character", default = NULL,
              help = "local MMCIF file (upload)", metavar = "character"),
  make_option(c("-w", "--workdir"), type = "character", default = NULL,
              help = "working directory", metavar = "character")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

if (is.null(opt$workdir)) {
  workDir = "."
} else {
  workDir = opt$workdir
}

setwd(workDir)

if(! is.null(opt$pdb)) {
  ## upload local PDB structure in '.pdb' format
  inputPDB = opt$pdb
} else if (! is.null(opt$id)) {
  ## download PDB structure based on PDB identifier
  get.pdb(opt$id, format = "cif", path = ".")
  pdbConversionName = sub("\\.cif$", ".pdb", opt$mmcif)
  command0 = paste("gemmi convert", opt$mmcif, pdbConversionName)
  system_status0 = system(command0)
  inputPDB = pdbConversionName
} else if (! is.null(opt$mmcif)) {
  ## convert from '.mmcif' format to '.pdb' format
  pdbConversionName = sub("\\.cif$", ".pdb", opt$mmcif)
  command0 = paste("gemmi convert", opt$mmcif, pdbConversionName)
  system_status0 = system(command0)
  inputPDB = pdbConversionName
} else {
  stop("No valid input. Get help with 'Rscript popscomp_standalone.R --help'.")
}


#________________________________________________________________________________
## ISO: split input PDB into chains
## 'pdbsplit' is a function from the 'bio3d' library
chain.files = pdbsplit(pdb.files = paste(inputPDB, sep = "/"),  path = ".", multi = FALSE)

## if input PDB is not a complex, return without any computations
##   bacause this routine "popscompR" is intended for processing protein complexes
if (length(chain.files) < 2) {
  stop("Single-chain structure: POPScomp not applicable")
}

chain.files.short = sub('\\.pdb$', '', basename(as.character(chain.files)))

#________________________________________________________________________________
## ISO: run POPS over all single (= isolated) chains via system (= shell) call
message("Isolated chains")
exit_codes = sapply(1:length(chain.files), function(x) {
 	command1 = paste0("pops --outDirName ", ".",
          " --rout --routPrefix ", paste0(chain.files.short[x], ".iso"),
					" --residueOut",
					" --pdb ", chain.files[x], " 1> ", chain.files.short[x], ".o",
					" 2> ", chain.files.short[x], ".e")
	print(command1)
	system_status1 = system(command1, wait = TRUE)
	message("  chain ", x, ": ", chain.files[x], "  exit code: ", system_status1)
})

## Concatenate output files of single (ISO = isolated) chains.
## Residue
command3 = paste0("tail -q -n+2 ", "*.iso.rpopsResidue >> ", "isoSASA.rpopsResidue")
system_status3 = system(command3, wait = TRUE)

#________________________________________________________________________________
## PAIR: create PDB files for all pairwise chain combinations
message("Paired chains")
pair.cmbn = combn(length(chain.files), 2)
chainpair.files = vector()
chainpair.files = sapply(1:dim(pair.cmbn)[2], function(x) {
 	## name of paired chain PDB file to create
 	chainpair.files[[x]] = paste0(chain.files.short[pair.cmbn[1, x]], "-",
                	              chain.files.short[pair.cmbn[2, x]], ".pdb")
 	## concatenate single chain PDB files to paired chain PDB files
	 command5 = paste("cat", chain.files[pair.cmbn[1, x]],
    	                     chain.files[pair.cmbn[2, x]], ">",
        	                 chainpair.files[[x]])
 	system_status5 = system(command5, wait = TRUE)
 	paste("  chain pair:", x, " exit code:", system_status5)
 	return(chainpair.files[[x]])
})

chainpair.files.short = sub('\\.pdb$', '', basename(as.character(chainpair.files)))

#________________________________________________________________________________
## PAIR: run POPS over all pairwise chain combinations via system (= shell) call
exit_codes = sapply(1:length(chainpair.files), function(x) {
	command6 = paste0("pops",
					" --rout --routPrefix ", paste0(chainpair.files.short[x], ".pair"),
					" --residueOut",
					" --pdb ", chainpair.files[x], " 1> ", "POPScomp_chainpair", x, ".o",
                    " 2> ", "POPScomp_chainpair", x, ".e")
	system_status6 = system(command6, wait = TRUE)
	message("  chain ", x, ": ", chainpair.files[x], "  exit code: ", system_status6)
})

#________________________________________________________________________________
## read SASA files
message("Processing SASA files")
## the data structure will be a list (levels = 'rpopsLevel') of lists (structures)
rpopsLevel = c("rpopsResidue")

message("SASA files of isolated chains")
## ISO: initialise list of lists with predefined number of output files
iso.sasa.level.files = vector(mode = "list", length = length(rpopsLevel))
iso.veclist = function(x) { vector(mode = "list", length = length(chain.files)) }
iso.sasa.level.files = lapply(iso.sasa.level.files, iso.veclist)

## read ISO SASA files
for (i in 1:length(chain.files)) {
	## read isolated chain output
	iso.sasa.level.files[[1]][[i]] = read.table(paste0(chain.files.short[i],
        	                                ".iso.", rpopsLevel[1]),
            	                            header = TRUE, stringsAsFactors = FALSE)
}
names(iso.sasa.level.files[[1]]) = chain.files.short

## PAIR: initialise list of lists with predefined number of output files
message("SASA files of paired chains")
pair.sasa.level.files = vector(mode = "list", length = length(rpopsLevel))
pair.veclist = function(x) { vector(mode = "list", length = dim(pair.cmbn)[2]) }
pair.sasa.level.files = lapply(pair.sasa.level.files, pair.veclist)

## read PAIR SASA files
for (i in 1:dim(pair.cmbn)[2]) {
	## read paired chain output
	pair.sasa.level.files[[1]][[i]] = read.table(paste0(chainpair.files.short[i],
      	                                ".pair.", rpopsLevel[1]),
          	                            header = TRUE, stringsAsFactors = FALSE)
}
names(pair.sasa.level.files[[1]]) = chainpair.files.short

#________________________________________________________________________________
## DIFF: compute SASA differences (POPScomp values)
## 'pair.cmbn' contains the order of PAIR files as column order and
##   the index of ISO files as column elements. That way the match between
##   PAIR and ISO files is reconstructed here.
## initialise list of lists with predefined number of SASA difference tables
#message("Compuing SASA differences")
diff.sasa.level = vector(mode = "list", length = length(rpopsLevel))
diff.veclist = function(x) { vector(mode = "list", length = dim(pair.cmbn)[2]) }
diff.sasa.level = lapply(diff.sasa.level, diff.veclist)

## compute SASA differences
iso.rbind.tmp = rbind(iso.sasa.level.files[[1]][[1]],
                      iso.sasa.level.files[[1]][[2]])

## using the log-ratio in a robust normalised form (Zrobust_Q)
#D_Phob.A.2 = round(iso.rbind.tmp[ , "Phob.A.2"] - pair.sasa.level.files[[1]][[1]][ , "Phob.A.2"], digits = 2)
#D_Phil.A.2 = round(iso.rbind.tmp[ , "Phil.A.2"] - pair.sasa.level.files[[1]][[1]][ , "Phil.A.2"], digits = 2)
#D_SASA.A.2 = round(iso.rbind.tmp[ , "SASA.A.2"] - pair.sasa.level.files[[1]][[1]][ , "SASA.A.2"], digits = 2)
logratio_Q.SASA = round(log2(iso.rbind.tmp[ , "Q.SASA."] / pair.sasa.level.files[[1]][[1]][ , "Q.SASA."]), digits = 2)
is.lix = logratio_Q.SASA > 0
#Z_Q.SASA = logratio_Q.SASA / sd(logratio_Q.SASA)
is.logratio_Q.SASA = logratio_Q.SASA[is.lix]
Zrobust_Q.SASA = (logratio_Q.SASA[is.lix] - median(logratio_Q.SASA[is.lix])) / (1.4862 * mad(logratio_Q.SASA[is.lix]))

## final residue selection is SASA <= 0.25 and Zrobust > 1
sel.lix = (pair.sasa.level.files[[1]][[1]][ , "Q.SASA."][is.lix] <= 0.25) & (Zrobust_Q.SASA > 1)
final.ix = (which(is.lix))[sel.lix]

#________________________________________________________________________________
## write results
write.table(pair.sasa.level.files[[1]][[1]][final.ix, ], file = "Qinterface.dat")


#===============================================================================
