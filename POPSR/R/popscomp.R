#! /usr/bin/R

#===============================================================================
# POPSR package
# Process complex structures
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
	inputPDB = "1f3r.pdb";
	outdir = tempdir();
	chain.files = pdbsplit(inputPDB,  path = outdir, multi = FALSE);
	## short names
	chain.files.short = sub('\\.pdb$', '', basename(chain.files));

	#________________________________________________________________________________
	## run POPS over all single chains via system (= shell) call
	sapply(1:length(chainpair.files), function(x) {
	  command = paste0("pops --outDirName ", outdir,
	                   " --rout --atomOut --residueOut --chainOut ",
	                   "--pdb ", chainpair.files[x], " 1> ", outdir, "/POPScomp_chainpair", x, ".o",
	                   " 2> ", outdir, "/POPScomp_chainpair", x, ".e");
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
}

#===============================================================================
