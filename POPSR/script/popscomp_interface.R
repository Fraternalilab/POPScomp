#!/usr/bin/env Rscript

#===============================================================================
# POPSR package
# popscomp_interface.R: Implementation of the POPSCOMP functionality,
# i.e. processing of complex structures to compute SASA difference values,
# and selection of the interface residues of each chain pair.
# Writes the POPS output files of the single chains and chain pairs
#   plus one '<chain pair>_Qinterface.dat' table per chain pair.
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

## input files are named relative to the directory the user started in,
## so they are resolved before the working directory is changed
if (! is.null(opt$pdb)) {
  opt$pdb = normalizePath(opt$pdb, mustWork = TRUE)
}
if (! is.null(opt$mmcif)) {
  opt$mmcif = normalizePath(opt$mmcif, mustWork = TRUE)
}

setwd(workDir)

if(! is.null(opt$pdb)) {
  ## upload local PDB structure in '.pdb' format
  ## a '.pdb' input needs no conversion; it is copied into the working directory
  ## so that the chain files and all output land in --workdir, not beside the source
  pdbConversionName = basename(opt$pdb)
  if (! identical(normalizePath(dirname(opt$pdb)), normalizePath("."))) {
    if (! file.copy(opt$pdb, pdbConversionName, overwrite = TRUE)) {
      stop("Could not copy ", opt$pdb, " into the working directory")
    }
  }
  inputPDB = pdbConversionName
} else if (! is.null(opt$id)) {
  ## download PDB structure based on PDB identifier
  if (! grepl("^[0-9A-Za-z]{4}$", opt$id)) {
    stop("Not a PDB identifier: ", opt$id)
  }
  get.pdb(opt$id, format = "cif", path = ".")
  ## 'get.pdb' writes '<id>.cif' into the working directory
  cifName = paste0(opt$id, ".cif")
  pdbConversionName = paste0(opt$id, ".pdb")
  command0 = paste("gemmi convert", shQuote(cifName), shQuote(pdbConversionName))
  system_status0 = system(command0)
  if (system_status0 != 0) {
    stop("Conversion of ", cifName, " failed with exit code ", system_status0)
  }
  inputPDB = pdbConversionName
} else if (! is.null(opt$mmcif)) {
  ## convert from '.mmcif' format to '.pdb' format;
  ## both '.cif' and '.cif.gz' are accepted and the input file is never written to
  pdbConversionName = sub("\\.cif(\\.gz)?$", ".pdb", basename(opt$mmcif))
  if (pdbConversionName == basename(opt$mmcif)) {
    stop("MMCIF input is expected to end in '.cif' or '.cif.gz': ", opt$mmcif)
  }
  if (grepl("\\.gz$", opt$mmcif)) {
    command0 = paste("zcat", shQuote(opt$mmcif), "|", "gemmi convert -", shQuote(pdbConversionName))
  } else {
    command0 = paste("gemmi convert", shQuote(opt$mmcif), shQuote(pdbConversionName))
  }
  system_status0 = system(command0)
  if (system_status0 != 0) {
    stop("Conversion of ", opt$mmcif, " failed with exit code ", system_status0)
  }
  inputPDB = pdbConversionName
} else {
  stop("No valid input. Get help with 'Rscript popscomp_interface.R --help'.")
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
          " --rout --routPrefix ", shQuote(paste0(chain.files.short[x], ".iso")),
					" --residueOut",
					" --pdb ", shQuote(chain.files[x]),
					" 1> ", shQuote(paste0(chain.files.short[x], ".o")),
					" 2> ", shQuote(paste0(chain.files.short[x], ".e")))
	system_status1 = system(command1, wait = TRUE)
	message("  chain ", x, ": ", chain.files[x], "  exit code: ", system_status1)
	return(system_status1)
})
if (any(exit_codes != 0)) {
	stop("POPS failed on isolated chain(s): ",
		paste(chain.files.short[exit_codes != 0], collapse = ", "))
}

## Concatenate output files of single (ISO = isolated) chains:
## the header line of the first chain, then the data lines of all chains.
## The file is truncated ('>'), so a re-run in the same directory does not
## append to the table of the previous run. Only the chains of this structure
## are listed, so left-over '*.iso.*' files of other structures are not swept in.
iso.residue.files = paste0(chain.files.short, ".iso.rpopsResidue")
command3 = paste("head -1", shQuote(iso.residue.files[1]), "> isoSASA.rpopsResidue &&",
                 "tail -q -n+2", paste(shQuote(iso.residue.files), collapse = " "),
                 ">> isoSASA.rpopsResidue")
system_status3 = system(command3, wait = TRUE)
if (system_status3 != 0) {
	stop("Concatenation of isolated-chain SASA files failed with exit code ", system_status3)
}

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
	 command5 = paste("cat", shQuote(chain.files[pair.cmbn[1, x]]),
    	                     shQuote(chain.files[pair.cmbn[2, x]]), ">",
        	                 shQuote(chainpair.files[[x]]))
 	system_status5 = system(command5, wait = TRUE)
 	paste("  chain pair:", x, " exit code:", system_status5)
 	return(chainpair.files[[x]])
})

chainpair.files.short = sub('\\.pdb$', '', basename(as.character(chainpair.files)))

#________________________________________________________________________________
## PAIR: run POPS over all pairwise chain combinations via system (= shell) call
exit_codes = sapply(1:length(chainpair.files), function(x) {
	command6 = paste0("pops --outDirName ", ".",
					" --rout --routPrefix ", shQuote(paste0(chainpair.files.short[x], ".pair")),
					" --residueOut",
					" --distMatCAOut ", shQuote(paste0(chainpair.files.short[x], ".distMatCA.out")),
					" --pdb ", shQuote(chainpair.files[x]),
					" 1> ", shQuote(paste0("POPScomp_chainpair", x, ".o")),
                    " 2> ", shQuote(paste0("POPScomp_chainpair", x, ".e")))
	system_status6 = system(command6, wait = TRUE)
	message("  chain pair ", x, ": ", chainpair.files[x], "  exit code: ", system_status6)
	return(system_status6)
})
if (any(exit_codes != 0)) {
	stop("POPS failed on chain pair(s): ",
		paste(chainpair.files.short[exit_codes != 0], collapse = ", "))
}

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
## DIFF: compute interface residues per chain pair
## 'pair.cmbn' contains the order of PAIR files as column order and
##   the index of ISO files as column elements. That way the match between
##   PAIR and ISO files is reconstructed here.
message("Computing interface residues")

for (i in 1:dim(pair.cmbn)[2]) {
	## the two isolated chains of this pair, in the order in which they were
	## concatenated into the pair structure
	iso.rbind.tmp = rbind(iso.sasa.level.files[[1]][[pair.cmbn[1, i]]],
	                      iso.sasa.level.files[[1]][[pair.cmbn[2, i]]])
	pair.tmp = pair.sasa.level.files[[1]][[i]]

	## The two tables are compared row by row, so they must describe the same
	## residues in the same order. That is verified here instead of assumed:
	## a silent mismatch would subtract unrelated residues from each other.
	res.key = function(x) paste(x[ , "Chain"], x[ , "ResidNr"], x[ , "iCode"], sep = ":")
	if (! identical(res.key(iso.rbind.tmp), res.key(pair.tmp))) {
		stop("Residues of the isolated chains and of the chain pair ",
			chainpair.files.short[i], " do not match")
	}

	## using the log-ratio in a robust normalised form (Zrobust_Q)
	## Residues that are fully buried in either structure have Q.SASA = 0, which
	## would give an infinite or undefined log-ratio; they are excluded here, as
	## are residues that lose no accessibility (log-ratio <= 0).
	q.iso = iso.rbind.tmp[ , "Q.SASA."]
	q.pair = pair.tmp[ , "Q.SASA."]
	logratio_Q.SASA = rep(NA_real_, length(q.iso))
	valid = is.finite(q.iso) & is.finite(q.pair) & (q.iso > 0) & (q.pair > 0)
	logratio_Q.SASA[valid] = round(log2(q.iso[valid] / q.pair[valid]), digits = 2)
	is.lix = ! is.na(logratio_Q.SASA) & (logratio_Q.SASA > 0)

	if (! any(is.lix)) {
		message("  chain pair ", chainpair.files.short[i], ": no buried residues")
		next
	}

	## robust z-score: 'mad' already scales by 1.4826 for consistency with the
	## standard deviation of a normal distribution, so it must not be scaled again
	lr = logratio_Q.SASA[is.lix]
	lr.mad = mad(lr)
	if (lr.mad == 0) {
		message("  chain pair ", chainpair.files.short[i],
			": log-ratios have zero deviation, no residue selected")
		next
	}
	Zrobust_Q.SASA = (lr - median(lr)) / lr.mad

	## final residue selection is Q.SASA <= 0.25 and Zrobust > 1
	sel.lix = (q.pair[is.lix] <= 0.25) & (Zrobust_Q.SASA > 1)
	final.ix = (which(is.lix))[sel.lix]

	#________________________________________________________________________________
	## write results: one table per chain pair
	outName = paste0(chainpair.files.short[i], "_Qinterface.dat")
	write.table(pair.tmp[final.ix, ], file = outName)
	message("  chain pair ", chainpair.files.short[i], ": ",
		length(final.ix), " interface residues -> ", outName)
}


#===============================================================================
