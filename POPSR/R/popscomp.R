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

library("parallel");

get.pdb("1f3r", path = "/tmp");
pdb = read.pdb("/tmp/1f3r.pdb");
setwd("/home/jkleinj/2COPIES/develop/POPSsuite/POPScomp/POPSR/inst/popsr");
outDir = ".";

#_______________________________________________________________________________
#' sasa : An S4 class for SASA data
#' @slot valueAtomPairchain: SASA values of paired chain at atom resolution
#' @slot valueAtomChain1: same for the first isolated chain
#' @slot valueAtomChain2: same for the second isolated chain
#' ... equally for the other resolutions 'Residue', 'Chain' and 'Molecule'
#' @slot diffAtomChain1: SASA differences of chain 1 at atom resolution
#' @slot diffAtomChain2: SASA differences of chain 2 at atom resolution
#' ... equally for the other resolutions 'Residue' and 'Chain'
sasa <- setClass(
  "sasa",

  slots = c(
    valueAtomPairchain = "data.frame",
    valueAtomChain1 = "data.frame",
    valueAtomChain2 = "data.frame",
    valueResiduePairchain = "data.frame",
    valueResidueChain1 = "data.frame",
    valueResidueChain2 = "data.frame",
    valueChainPairchain = "data.frame",
    valueChainChain1 = "data.frame",
    valueChainChain2 = "data.frame",
    valueMoleculePairchain = "data.frame",
    valueMoleculeChain1 = "data.frame",
    valueMoleculeChain2 = "data.frame",
    diffAtomChain1 = "data.frame",
    diffAtomChain2 = "data.frame",
    diffResidueChain1 = "data.frame",
    diffResidueChain2 = "data.frame",
    diffChainChain1 = "data.frame",
    diffChainChain2 = "data.frame"
  )
)

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
	nCore = detectCores() - 1;

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
	## first just the header line
	command0 = paste0("head -n 1 ", outDir, "/id.rpopsAtom > ", outDir, "/isoSASA.rpopsAtom");
	system_status0 = system(command0);
	## then concatenate all putput files, without file name header (-q) and table header (-n+2)
	command1 = paste0("tail -q -n+2 ", outDir, "/*.iso.rpopsAtom >> ", outDir, "/isoSASA.rpopsAtom");
	system_status1 = system(command1);
	command2 = paste0("tail -q -n+2 ", outDir, "/*.iso.rpopsResidue >> ", outDir, "/isoSASA.rpopsResidue");
	system_status2 = system(command2);
	command3 = paste0("tail -q -n+2 ", outDir, "/*.iso.rpopsChain >> ", outDir, "/isoSASA.rpopsChain");
	system_status3 = system(command3);

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
  ## initialise list of lists with predefined number of output files
  diff.sasa.level = vector(mode = "list", length = length(rpopsLevel));
  diff.veclist = function(x) { vector(mode = "list", length = length(dim(pair.cmbn)[2])) };
  diff.sasa.level = lapply(diff.sasa.level, diff.veclist);
  sasadiff.level.tables = lapply(1:(length(rpopsLevel) - 1), function(y) {
  ## compute SASA differences
  for (j in 1:length(rpopsLevel)) {
    for (i in 1:length(dim(pair.cmbn)[2])) {
	    ## rbind ISO chain SASAs
      iso.tmp = rbind(iso.sasa.level.files[[j]][[pair.cmbn[1, i]]],
                      iso.sasa.level.files[[j]][[pair.cmbn[2, i]]]);
      ## assert consistency between 'rbind' ISO files and PAIR file
      stopifnot(dim(iso.tmp) == dim(pair.sasa.level.files[[j]][[i]]));
      ## DIFF values
      diff.tmp = iso.tmp[ , "SASA.A.2"] - pair.sasa.level.files[[j]][[i]][ , "SASA.A.2"];
	  });
	});

	return(list(sasapair.level.files, sasadiff.level.tables));
}

#===============================================================================
