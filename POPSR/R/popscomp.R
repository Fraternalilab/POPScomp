#! /usr/bin/R

#===============================================================================
# POPSR package
# popscompR: Processing of complex structures to compute SASA difference values
# Returns a list of POPS output files of single-chain and pair-chain structures
#   plus a list of per-chain difference values.
# (C) 2019 Jens Kleinjung and Franca Fraternali
#===============================================================================

library("bio3d");
## for atom selection mechanisms see:
## http://thegrantlab.org/bio3d/tutorials/structure-analysis

library("parallel");

popscompR = function(inputPDB, outdir) {
	## number of cores
	nCore = detectCores() - 1;

	#________________________________________________________________________________
	## split input PDB into chains
	chain.files = pdbsplit(paste(outdir, inputPDB, sep = "/"),  path = outdir, multi = FALSE);
	## short names
	chain.files.short = sub('\\.pdb$', '', basename(chain.files));

	#________________________________________________________________________________
	## run POPS over all single chains via system (= shell) call
	sapply(1:length(chain.files), function(x) {
	  command = paste0("pops --outDirName ", outdir,
	                   " --rout --atomOut --residueOut --chainOut ",
	                   "--pdb ", chain.files[x], " 1> ", outdir, "/POPScomp_chain", x, ".o",
	                   " 2> ", outdir, "/POPScomp_chain", x, ".e");
	  system_status = system(command);
	  paste("Exit code:", system_status);
	});

	#________________________________________________________________________________
	## create PDB files for all pairwise chain combinations
	pair.cmbn = combn(length(chain.files), 2);
	chainpair.files = list();
	chainpair.files = sapply(1:dim(pair.cmbn)[2], function(x) {
	  chainpair.files[[x]] = paste0(outdir, "/",
	                          chain.files.short[pair.cmbn[1, x]], "-",
	                          chain.files.short[pair.cmbn[2, x]], ".pdb");
	  command = paste("cat", chain.files[pair.cmbn[1, x]],
	                         chain.files[pair.cmbn[2, x]], ">",
	                         chainpair.files[[x]]);
	  system_status = system(command);
	  paste("Chain pair:", x, "; Exit code:", system_status);
	  return(chainpair.files);
	});

	#________________________________________________________________________________
	## run POPS over all pairwise chain combinations via system (= shell) call
	sapply(1:length(chainpair.files), function(x) {
	  command = paste0("pops --outDirName ", outdir,
	                  " --rout --atomOut --residueOut --chainOut ",
	                  "--pdb ", chainpair.files[x], " 1> ", outdir, "/POPScomp_chainpair", x, ".o",
                                                  " 2> ", outdir, "/POPScomp_chainpair", x, ".e");
	  system_status = system(command);
	  paste("Exit code:", system_status);
	});

	#________________________________________________________________________________
	## read SASA files
  ## the structure will be a list (all chain pairs) of lists (chainpair12, chain1, chain2)
	rpopsLevel = c(".rpopsAtom", ".rpopsResidue", ".rpopsChain", ".rpopsMolecule");
  ## sasapair.level is a list with elements at resolutions: atom, residue, chain, molecule
	sasapair.level.files = lapply(1:length(rpopsLevel), function(y) {
	  ## 'sasapair' is a list with a chain pair SASA and the two single chain SASAs as elements
	  sasapair.files = lapply(1:dim(pair.cmbn)[2], function(x) {
      sasa.files = list();
      sasa.files[[1]] = read.table(paste0(chainpair.files[x], rpopsLevel[y]), header = TRUE);
      sasa.files[[2]] = read.table(paste0(outdir, "/", chain.files.short[pair.cmbn[1, x]],
                                          ".pdb", rpopsLevel[y]), header = TRUE);
      sasa.files[[3]] = read.table(paste0(outdir, "/", chain.files.short[pair.cmbn[2, x]],
                                          ".pdb", rpopsLevel[y]), header = TRUE);
      names(sasa.files) = c(paste0(sasa.files[[2]][1, "Chain"], sasa.files[[3]][1, "Chain"]),
                            as.character(sasa.files[[2]][1, "Chain"]),
                            as.character(sasa.files[[3]][1, "Chain"]));
      return(sasa.files);
    });
  });
	names(sasapair.level.files) = rpopsLevel;

	#________________________________________________________________________________
	## compute SASA differences (POPScomp values)
	## merging first the single-chain with the pair-chain values
	## by using the 'AtomNr' (atom serial number), because that refers to the line of
	##   the coordinate entry and is therefore unique and contiguous
	bynameLevel = c("AtomNr", "Chain", "Id");
	sasapair.level.files = lapply(1:(length(rpopsLevel) - 1), function(y) {
	  ## sasadiff is a list of SASA difference values for two chains as list elements
		sasadiff.tables = lapply(1:dim(pair.cmbn)[2], function(x) {
	    sasa.diff = list();
	    ## SASA difference first chain
	    sasa.diff[[1]] = merge(sasapair.level.files[[y]][[x]][[2]], sasapair.level.files[[y]][[x]][[1]],
	                           by = bynameLevel[y], all = FALSE);
	    sasa.diff[[1]]$diff = sasa.diff[[1]]$SASA.A.2.x - sasa.diff[[1]]$SASA.A.2.y;
	    ## SASA difference second chain
	    sasa.diff[[2]] = merge(sasapair.level.files[[y]][[x]][[3]], sasapair.level.files[[y]][[x]][[1]],
	                           by = bynameLevel[y], all = FALSE);
	    sasa.diff[[2]]$diff = sasa.diff[[2]]$SASA.A.2.x - sasa.diff[[2]]$SASA.A.2.y;
	  });
	});

	return(list(sasapair.files[[1]], sasadiff.tables[[1]]));
}

#===============================================================================
