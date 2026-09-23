#! /usr/bin/R

#===============================================================================
# POPSR package
# popscomp.R: Implementation of the POPSCOMP functionality,
# i.e. processing of complex structures to compute SASA difference values.
# Returns a list of POPS output files for single-chain and pair-chain structures
#   plus a list of buried SASA values.
#
# (C) 2019-2026 Jens Kleinjung and Franca Fraternali
#===============================================================================

## for atom selection mechanisms see:
## http://thegrantlab.org/bio3d/tutorials/structure-analysis

#_______________________________________________________________________________
## POPScomp function implemented in R
## The following prefixes are used to label the sections and output files
##   for clarity (DIFF output files are called 'delta' for historic reasons):
## ID: the default '--popsr' prefix of POPS for the unmodified input PDB
##      (computed by 'input$popscomp' function in 'app.R')
## ISO: POPS on isolated chains
## PAIR: POPS on paired chains
## DIFF: difference between sum of isolated chain SASA and paired chain SASA
popscompR = function(inputPDB, outDir, coarse = FALSE) {

	## Path of the POPS program: the POPS_BIN environment variable, else the
	## program on the PATH, else the location in the Shiny Docker image.
	pops_bin = Sys.getenv("POPS_BIN", unset = "")
	if (! nzchar(pops_bin)) {
		pops_bin = Sys.which("pops")
	}
	if (! nzchar(pops_bin)) {
		pops_bin = "/build/install/usr/local/bin/pops"
	}
	if (! file.exists(pops_bin)) {
		stop("POPS program not found: set POPS_BIN or put 'pops' on the PATH");
	}

	## atomistic or coarse-grained (one centre per residue) computation
	popsMode = if (coarse) " --coarse" else "";
	## the atom level exists only in atomistic mode
	rpopsLevel = c("rpopsAtom", "rpopsResidue", "rpopsChain", "rpopsMolecule");
	if (coarse) {
		rpopsLevel = c("rpopsResidue", "rpopsChain", "rpopsMolecule");
	}

	#________________________________________________________________________________
	## ISO: split input PDB into chains
	chain.files = pdbsplit(paste(outDir, inputPDB, sep = "/"),  path = outDir, multi = FALSE);

	## if input PDB is not a complex, return without any computations
	##   bacause this routine "popscompR" is intended for processing protein complexes
	if (length(chain.files) < 2) {
	  message("Single-chain POPScomp");
	  return(0);
	}

	chain.files.short = sub('\\.pdb$', '', basename(as.character(chain.files)));

	#________________________________________________________________________________
	## ISO: run POPS over all single (= isolated) chains via system (= shell) call
	## File names are quoted: they come from the input structure and may contain
	##   spaces or characters that the shell would otherwise interpret.
	iso.ok = sapply(1:length(chain.files), function(x) {
	  command = paste0(shQuote(pops_bin), " --outDirName ", shQuote(outDir),
	                   " --rout --routPrefix ", shQuote(paste0(chain.files.short[x], ".iso")),
	                   " --atomOut --residueOut --chainOut", popsMode,
	                   " --pdb ", shQuote(chain.files[x]),
	                   " 1> ", shQuote(paste0(outDir, "/", chain.files.short[x], ".o")),
	                   " 2> ", shQuote(paste0(outDir, "/", chain.files.short[x], ".e")));
	  system_status = system(command, wait = TRUE);
	  ## POPS declines structures whose surface is not well defined (a chain of
	  ##   water or ligand residues, for example). Such a chain is dropped with a
	  ##   warning, so that the interfaces of the remaining chains are still computed.
	  ok = (system_status == 0) &&
	       all(file.exists(paste0(outDir, "/", chain.files.short[x], ".iso.", rpopsLevel)));
	  if (! ok) {
	    warning("POPS could not process chain ", chain.files.short[x],
	            " (exit code ", system_status, "); the chain is skipped", call. = FALSE);
	  }
	  return(ok);
	});

	chain.files = chain.files[iso.ok];
	chain.files.short = chain.files.short[iso.ok];
	if (length(chain.files) < 2) {
	  message("Fewer than two usable chains: POPScomp not applicable");
	  return(0);
	}

	## Concatenate output files of single (ISO = isolated) chains.
	## We do that here because there is only one tab on the interface for each resolution level
	##   and the number of chains to be processed/shown will vary. Otherwise we would need
	##   a dynamic tab structure on the interface that creates a tab for each chain.
	## The header is taken from the first chain file of the same level, so that the
	##   column names always match the table below them, and the files of this run are
	##   listed explicitly, so that neither the shell's alphabetical glob order nor
	##   left-over files of an earlier run can enter the table.
	for (lvl in setdiff(rpopsLevel, "rpopsMolecule")) {
	  iso.paths = paste0(outDir, "/", chain.files.short, ".iso.", lvl);
	  out.path = paste0(outDir, "/isoSASA.", lvl);
	  header = readLines(iso.paths[1], n = 1);
	  rows = unlist(lapply(iso.paths, function(f) readLines(f)[-1]));
	  writeLines(c(header, rows), out.path);
	}

	#________________________________________________________________________________
	## PAIR: create PDB files for all pairwise chain combinations
	pair.cmbn = combn(length(chain.files), 2);
	chainpair.files = sapply(1:dim(pair.cmbn)[2], function(x) {
	  ## name of paired chain PDB file to create
	  pair.file = paste0(outDir, "/",
	                     chain.files.short[pair.cmbn[1, x]], "-",
	                     chain.files.short[pair.cmbn[2, x]], ".pdb");
	  ## concatenate single chain PDB files to paired chain PDB files
	  writeLines(c(readLines(chain.files[pair.cmbn[1, x]]),
	               readLines(chain.files[pair.cmbn[2, x]])), pair.file);
	  return(pair.file);
	});

	chainpair.files.short = sub('\\.pdb$', '', basename(as.character(chainpair.files)));

	#________________________________________________________________________________
	## PAIR: run POPS over all pairwise chain combinations via system (= shell) call
	pair.ok = sapply(1:length(chainpair.files), function(x) {
	  command = paste0(shQuote(pops_bin), " --outDirName ", shQuote(outDir),
	                  " --rout --routPrefix ", shQuote(paste0(chainpair.files.short[x], ".pair")),
	                  " --atomOut --residueOut --chainOut", popsMode,
	                  " --pdb ", shQuote(chainpair.files[x]),
	                  " 1> ", shQuote(paste0(outDir, "/POPScomp_chainpair", x, ".o")),
	                  " 2> ", shQuote(paste0(outDir, "/POPScomp_chainpair", x, ".e")));
	  system_status = system(command, wait = TRUE);
	  ok = (system_status == 0) &&
	       all(file.exists(paste0(outDir, "/", chainpair.files.short[x], ".pair.", rpopsLevel)));
	  if (! ok) {
	    warning("POPS could not process chain pair ", chainpair.files.short[x],
	            " (exit code ", system_status, "); the pair is skipped", call. = FALSE);
	  }
	  return(ok);
	});

	pair.cmbn = pair.cmbn[ , pair.ok, drop = FALSE];
	chainpair.files.short = chainpair.files.short[pair.ok];
	if (dim(pair.cmbn)[2] == 0) {
	  message("No usable chain pair: POPScomp not applicable");
	  return(0);
	}

	#________________________________________________________________________________
	## read SASA files
	## the data structure will be a list (levels = 'rpopsLevel') of lists (structures)
	## Columns that identify a residue, atom or chain are read as character:
	##   'T' and 'F' are valid chain identifiers and would otherwise be read as
	##   the logical values TRUE and FALSE.
	## The column classes are set before reading, not corrected afterwards:
	##   read.table would already have turned a chain identifier 'T' or 'F' into
	##   the logical values TRUE and FALSE, which cannot be undone.
	read.rpops = function(path) {
	  header = scan(path, what = "", nlines = 1, quiet = TRUE);
	  colClasses = rep(NA_character_, length(header));
	  colClasses[header %in% c("AtomNe", "ResidNe", "Chain", "Id", "iCode",
	                           "AtomRange", "ResidRange")] = "character";
	  read.table(path, header = TRUE, stringsAsFactors = FALSE, colClasses = colClasses);
	}

	## ISO: initialise list of lists with predefined number of output files
	iso.sasa.level.files = vector(mode = "list", length = length(rpopsLevel));
	iso.veclist = function(x) { vector(mode = "list", length = length(chain.files)) };
	iso.sasa.level.files = lapply(iso.sasa.level.files, iso.veclist);

	## read ISO SASA files
	for (j in 1:length(rpopsLevel)) {
	  for (i in 1:length(chain.files)) {
	    ## read isolated chain output
	    iso.sasa.level.files[[j]][[i]] = read.rpops(paste0(outDir, "/", chain.files.short[i],
	                                        ".iso.", rpopsLevel[j]));
	  };
	  names(iso.sasa.level.files[[j]]) = chain.files.short;
	};
	names(iso.sasa.level.files) = rpopsLevel;

	## PAIR: initialise list of lists with predefined number of output files
	pair.sasa.level.files = vector(mode = "list", length = length(rpopsLevel));
	pair.veclist = function(x) { vector(mode = "list", length = dim(pair.cmbn)[2]) };
	pair.sasa.level.files = lapply(pair.sasa.level.files, pair.veclist);

	## read PAIR SASA files
	for (j in 1:length(rpopsLevel)) {
	  for (i in 1:dim(pair.cmbn)[2]) {
	    ## read paired chain output
	    pair.sasa.level.files[[j]][[i]] = read.rpops(paste0(outDir, "/", chainpair.files.short[i],
	                                        ".pair.", rpopsLevel[j]));
	  }
	  names(pair.sasa.level.files[[j]]) = chainpair.files.short;
	}
	names(pair.sasa.level.files) = rpopsLevel;

	#________________________________________________________________________________
	## DIFF: compute SASA differences (POPScomp values)
	## 'pair.cmbn' contains the order of PAIR files as column order and
	##   the index of ISO files as column elements. That way the match between
	##   PAIR and ISO files is reconstructed here.
	## initialise list of lists with predefined number of SASA difference tables
	diff.sasa.level = vector(mode = "list", length = length(rpopsLevel));
	diff.veclist = function(x) { vector(mode = "list", length = dim(pair.cmbn)[2]) };
	diff.sasa.level = lapply(diff.sasa.level, diff.veclist);

	## The isolated-chain rows and the chain-pair rows are subtracted from each
	##   other row by row, so they must describe the same atoms, residues or chains
	##   in the same order. The key is compared instead of assumed: a mismatch
	##   would silently subtract unrelated entries.
	entry.key = function(x, level) {
	  cols = switch(level,
	    rpopsAtom = c("AtomNe", "ResidNe", "Chain", "ResidNr", "iCode"),
	    rpopsResidue = c("ResidNe", "Chain", "ResidNr", "iCode"),
	    rpopsChain = c("Id"),
	    NULL);
	  if (is.null(cols)) return(NULL);
	  do.call(paste, c(lapply(cols, function(cn) x[[cn]]), sep = ":"));
	}

	## compute SASA differences
	for (j in 1:length(rpopsLevel)) {
	  level = rpopsLevel[j];
	  for (i in 1:dim(pair.cmbn)[2]) {
	    iso.rbind.tmp = rbind(iso.sasa.level.files[[j]][[pair.cmbn[1, i]]],
	                          iso.sasa.level.files[[j]][[pair.cmbn[2, i]]]);
	    pair.tmp = pair.sasa.level.files[[j]][[i]];

	    ## not for the molecule level, which is a single row per structure
	    if (level != "rpopsMolecule") {
	      ## assert consistency between 'rbind' ISO files and PAIR file
	      stopifnot(dim(iso.rbind.tmp) == dim(pair.tmp));
	      if (! identical(entry.key(iso.rbind.tmp, level), entry.key(pair.tmp, level))) {
	        stop("Entries of the isolated chains and of the chain pair ",
	             chainpair.files.short[i], " do not match at level ", level);
	      }
	      ## SASA DIFF values, applies to all levels
	      D_SASA.A.2 = round(iso.rbind.tmp[ , "SASA.A.2"] - pair.tmp[ , "SASA.A.2"], 2);
	      ## the apolar/polar split exists at residue and chain level, not per atom
	      if (all(c("Phob.A.2", "Phil.A.2") %in% colnames(pair.tmp))) {
	        D_Phob.A.2 = round(iso.rbind.tmp[ , "Phob.A.2"] - pair.tmp[ , "Phob.A.2"], digits = 2);
	        D_Phil.A.2 = round(iso.rbind.tmp[ , "Phil.A.2"] - pair.tmp[ , "Phil.A.2"], digits = 2);
	        ## an entry is part of the interface if any of the three areas changed:
	        ##   a buried entry whose apolar and polar changes cancel has D_SASA = 0
	        buried = (D_SASA.A.2 > 0) | (D_Phob.A.2 != 0) | (D_Phil.A.2 != 0);
	      } else {
	        buried = (D_SASA.A.2 > 0);
	      }
	    }

	    ## more level-specific delta values
	    if (level == "rpopsAtom") {
	      ## the identifying columns are taken from the chain-pair table, whose
	      ##   chain index and ranges refer to the complex
	      diff.tmp.df = cbind(pair.tmp, D_SASA.A.2);
	      diff.sasa.level[[j]][[i]] = diff.tmp.df[buried,
	        c("AtomNr", "AtomNe", "ResidNe", "Chain", "ResidNr", "iCode",
	          "D_SASA.A.2", "AtomTp", "AtomGp")];
	    } else if (level == "rpopsResidue") {
	      diff.tmp.df = cbind(pair.tmp, D_Phob.A.2, D_Phil.A.2, D_SASA.A.2);
	      diff.sasa.level[[j]][[i]] = diff.tmp.df[buried,
	        c("ResidNe", "Chain", "ResidNr", "iCode", "D_Phob.A.2", "D_Phil.A.2", "D_SASA.A.2")];
	    } else if (level == "rpopsChain") {
	      diff.tmp.df = cbind(pair.tmp, D_Phob.A.2, D_Phil.A.2, D_SASA.A.2);
	      diff.sasa.level[[j]][[i]] = diff.tmp.df[buried,
	        c("Chain", "Id", "AtomRange", "ResidRange", "D_Phob.A.2", "D_Phil.A.2", "D_SASA.A.2")];
	    } else if (level == "rpopsMolecule") {
	      diff.tmp.df = iso.sasa.level.files[[j]][[pair.cmbn[1, i]]] +
	                    iso.sasa.level.files[[j]][[pair.cmbn[2, i]]] -
	                    pair.tmp;
	      diff.tmp.df = round(diff.tmp.df[ , c("Phob.A.2", "Phil.A.2", "SASA.A.2")], digits = 2);
	      ## chains that do not touch differ only by rounding noise
	      diff.tmp.df[abs(diff.tmp.df) < 0.05] = 0;
	      colnames(diff.tmp.df) = c("D_Phob.A.2", "D_Phil.A.2", "D_SASA.A.2");
	      diff.sasa.level[[j]][[i]] = diff.tmp.df;
	    }
	  };
	  names(diff.sasa.level[[j]]) = chainpair.files.short;
	};
	names(diff.sasa.level) = rpopsLevel;

	#________________________________________________________________________________
	## write DIFF SASA result files
	## ID SASA files have been created in the App
	## PAIR SASA files have been created here earlier
	## that completes the set of three types of output files
	for (j in 1:length(rpopsLevel)) {
	  diff.all = do.call(rbind, diff.sasa.level[[j]]);
	  if (is.null(diff.all)) {
	    next;
	  }
	  write.table(diff.all, paste0(outDir, "/", "deltaSASA.", rpopsLevel[j]));
	}

	return(diff.sasa.level);
}

#===============================================================================
