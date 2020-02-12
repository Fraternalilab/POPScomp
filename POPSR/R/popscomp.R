#! /usr/bin/R

#===============================================================================
# POPSR package
# popscomp.R: Implementation of the POPSCOMP functionality,
# i.e. processing of complex structures to compute SASA difference values.
# Returns a list of POPS output files for single-chain and pair-chain structures
#   plus a list of buried SASA values.
#
# (C) 2019-2020 Jens Kleinjung and Franca Fraternali
#===============================================================================

library("bio3d");
## for atom selection mechanisms see:
## http://thegrantlab.org/bio3d/tutorials/structure-analysis

#library("parallel");

#_______________________________________________________________________________
## POPScomp function implemented in R
## Prefixes:
## ID: the default '--popsr' prefix of POPS for the unmodified input PDB
##      (computed by 'input$popscomp' function in 'app.R')
## ISO: POPS on isolated chains
## PAIR: POPS on paired chains
## DIFF: difference between sum of isolated chain SASA and paired chain SASA
popscompR = function(inputPDB, outDir) {
	## number of cores
	#nCore = detectCores() - 1;

	#________________________________________________________________________________
	## ISO: split input PDB into chains
	chain.files = pdbsplit(paste(outDir, inputPDB, sep = "/"),  path = outDir, multi = FALSE);
	## short names
	chain.files.short = sub('\\.pdb$', '', basename(chain.files));

	#________________________________________________________________________________
	## ISO: run POPS over all single (= isolated) chains via system (= shell) call
	sapply(1:length(chain.files), function(x) {
	  command = paste0("pops --outDirName ", outDir,
	                   " --rout --routPrefix ", paste0(chain.files.short[x], ".iso"),
	                   " --atomOut --residueOut --chainOut",
	                   " --pdb ", chain.files[x], " 1> ", outDir, "/", chain.files.short[x], ".o",
	                   " 2> ", outDir, "/", chain.files.short[x], ".e");
	  system_status = system(command);
	  paste("Exit code:", system_status);
	});

	## Concatenate output files of single (= isolated) chains.
	## We do that here because there is only one tab on the interface for each resolution level
	##   and the number of chains to be processed/shown will vary. Otherwise we would need
	##   a dynamic tab structure on the interface that creates a tab for each chain.
	## Atom
	## first just the header line
	command1.1 = paste0("head -n 1 ", outDir, "/id.rpopsAtom > ", outDir, "/isoSASA.rpopsAtom");
	system_status1.1 = system(command1.1);
	## then concatenate output file, without file name header (-q) and table header (-n+2)
	command1.2 = paste0("tail -q -n+2 ", outDir, "/*.iso.rpopsAtom >> ", outDir, "/isoSASA.rpopsAtom");
	system_status1.2 = system(command1.2);
	## Residue
	command2.1 = paste0("head -n 1 ", outDir, "/id.rpopsResidue > ", outDir, "/isoSASA.rpopsResidue");
	system_status2.1 = system(command2.1);
	command2.2 = paste0("tail -q -n+2 ", outDir, "/*.iso.rpopsResidue >> ", outDir, "/isoSASA.rpopsResidue");
	system_status2.2 = system(command2.2);
	## Chain
	command3.1 = paste0("head -n 1 ", outDir, "/id.rpopsChain > ", outDir, "/isoSASA.rpopsChain");
	system_status3.1 = system(command3.1);
	command3.2 = paste0("tail -q -n+2 ", outDir, "/*.iso.rpopsChain >> ", outDir, "/isoSASA.rpopsChain");
	system_status3.2 = system(command3.2);

	#________________________________________________________________________________
	## PAIR: create PDB files for all pairwise chain combinations
	pair.cmbn = combn(length(chain.files), 2);
	chainpair.files = vector();
	chainpair.files = sapply(1:dim(pair.cmbn)[2], function(x) {
	  ## name of paired chain PDB file to create
	  chainpair.files[[x]] = paste0(outDir, "/",
	                          chain.files.short[pair.cmbn[1, x]], "-",
	                          chain.files.short[pair.cmbn[2, x]], ".pdb");
	  ## concatenate single chain PDB files to paired chain PDB files
	  command = paste("cat", chain.files[pair.cmbn[1, x]],
	                         chain.files[pair.cmbn[2, x]], ">",
	                         chainpair.files[[x]]);
	  system_status = system(command);
	  paste("Chain pair:", x, "; Exit code:", system_status);
	  return(chainpair.files);
	});

	chainpair.files.short = sub('\\.pdb$', '', basename(chainpair.files));

	#________________________________________________________________________________
	## PAIR: run POPS over all pairwise chain combinations via system (= shell) call
	sapply(1:length(chainpair.files), function(x) {
	  command = paste0("pops --outDirName ", outDir,
	                  " --rout --routPrefix ", paste0(chainpair.files.short[x], ".pair"),
	                  " --atomOut --residueOut --chainOut",
	                  " --pdb ", chainpair.files[x], " 1> ", outDir, "/POPScomp_chainpair", x, ".o",
                                                  " 2> ", outDir, "/POPScomp_chainpair", x, ".e");
	  system_status = system(command);
	  paste("Exit code:", system_status);
	});

	#________________________________________________________________________________
	## read SASA files
  ## the structure will be a list (levels) of lists (structures)
	rpopsLevel = c("rpopsAtom", "rpopsResidue", "rpopsChain", "rpopsMolecule");

	## initialise list of lists with predefined number of output files
	iso.sasa.level.files = vector(mode = "list", length = length(rpopsLevel));
	iso.veclist = function(x) { vector(mode = "list", length = length(chain.files)) };
	iso.sasa.level.files = lapply(iso.sasa.level.files, iso.veclist);
	## read ISO SASA files
	for (j in 1:length(rpopsLevel)) {
	  for (i in 1:length(chain.files)) {
	    ## read isolated chain output
	    iso.sasa.level.files[[j]][[i]] = read.table(paste0(outDir, "/", chain.files.short[i],
	                                        ".iso.", rpopsLevel[j]), header = TRUE);
	  };
	  names(iso.sasa.level.files[[j]]) = chain.files.short;
	};
  names(iso.sasa.level.files) = rpopsLevel;

  ## initialise list of lists with predefined number of output files
  pair.sasa.level.files = vector(mode = "list", length = length(rpopsLevel));
  pair.veclist = function(x) { vector(mode = "list", length = length(dim(pair.cmbn)[2])) };
  pair.sasa.level.files = lapply(pair.sasa.level.files, pair.veclist);
  ## read PAIR SASA files
  for (j in 1:length(rpopsLevel)) {
    for (i in 1:length(dim(pair.cmbn)[2])) {
      ## read paired chain output
      pair.sasa.level.files[[j]][[i]] = read.table(paste0(outDir, "/", chainpair.files.short[i],
                                          ".pair.", rpopsLevel[j]), header = TRUE);
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
  diff.veclist = function(x) { vector(mode = "list", length = length(dim(pair.cmbn)[2])) };
  diff.sasa.level = lapply(diff.sasa.level, diff.veclist);
  ## compute SASA differences
  for (j in 1:length(rpopsLevel)) {
    for (i in 1:length(dim(pair.cmbn)[2])) {
      if (j %in% 1:3) {
        ## rbind ISO chain SASAs
        iso.rbind.tmp = rbind(iso.sasa.level.files[[j]][[pair.cmbn[1, i]]],
                              iso.sasa.level.files[[j]][[pair.cmbn[2, i]]]);
        ## assert consistency between 'rbind' ISO files and PAIR file
        #stopifnot(dim(iso.rbind.tmp) == dim(pair.sasa.level.files[[j]][[i]]));
        print(paste(j, i, dim(iso.rbind.tmp), dim(pair.sasa.level.files[[j]][[i]])));
        ## SASA DIFF values, applies to all levels
        D_SASA.A.2 = round(iso.rbind.tmp[ , "SASA.A.2"] - pair.sasa.level.files[[j]][[i]][ , "SASA.A.2"], 2);
        ## more level-specific delta values
        if (j == 1) {
          diff.tmp.df = cbind(iso.rbind.tmp, D_SASA.A.2);
          diff.sasa.level[[j]][[i]] = diff.tmp.df[diff.tmp.df[ , "D_SASA.A.2"] > 0,
            c("ResidNe", "Chain", "ResidNr", "iCode", "D_SASA.A.2")];
        } else if (j == 2) {
          D_Phob.A.2 = round(iso.rbind.tmp[ , "Phob.A.2"] - pair.sasa.level.files[[j]][[i]][ , "Phob.A.2"], digits = 2);
          D_Phil.A.2 = round(iso.rbind.tmp[ , "Phil.A.2"] - pair.sasa.level.files[[j]][[i]][ , "Phil.A.2"], digits = 2);
          diff.tmp.df = cbind(iso.rbind.tmp, D_Phob.A.2, D_Phil.A.2, D_SASA.A.2);
          diff.sasa.level[[j]][[i]] = diff.tmp.df[diff.tmp.df[ , "D_SASA.A.2"] > 0,
            c("ResidNe", "Chain", "ResidNr", "iCode", "D_Phob.A.2", "D_Phil.A.2", "D_SASA.A.2")];
        } else if (j == 3) {
          D_Phob.A.2 = round(iso.rbind.tmp[ , "Phob.A.2"] - pair.sasa.level.files[[j]][[i]][ , "Phob.A.2"], digits = 2);
          D_Phil.A.2 = round(iso.rbind.tmp[ , "Phil.A.2"] - pair.sasa.level.files[[j]][[i]][ , "Phil.A.2"], digits = 2);
          diff.tmp.df = cbind(iso.rbind.tmp, D_Phob.A.2, D_Phil.A.2, D_SASA.A.2);
          diff.sasa.level[[j]][[i]] = diff.tmp.df[diff.tmp.df[ , "D_SASA.A.2"] > 0,
            c("Chain", "Id", "AtomRange", "ResidRange", "D_Phob.A.2", "D_Phil.A.2", "D_SASA.A.2")];
        }
      }
      if (j == 4) {
        diff.tmp.df = round(iso.sasa.level.files[[j]][[pair.cmbn[1, i]]] +
                            iso.sasa.level.files[[j]][[pair.cmbn[2, i]]] -
                            pair.sasa.level.files[[j]][[i]],
                            digits = 2);
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
    write.table(do.call(rbind, diff.sasa.level[[j]]), paste0(outDir, "/", "deltaSASA.", rpopsLevel[j]));
  }
}

#===============================================================================
