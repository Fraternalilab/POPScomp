#! /usr/bin/R

#===============================================================================
# POPSCOMP in R, first version specific for PKM2 
#===============================================================================

library("bio3d");
## for atom selection mechanisms see:
## http://thegrantlab.org/bio3d/tutorials/structure-analysis

library("parallel");
## number of cores
nCore = detectCores() - 1;

## all PDB structures
inFileDir = ".";
inFileName = list.files(path = inFileDir, full.names = FALSE, pattern = 'centroid.pdb$');

## general output directory
outDir = "popscomp";
dir.create(file.path(".", outDir), showWarnings = FALSE);

#_______________________________________________________________________________
## process domains 
for (i in 1:length(inFileName)) {
	pdb = read.pdb(paste(inFileDir, inFileName[i], sep = "/"));

	## structure-specific output directory
	outStrDir = substr(inFileName[i], 1, nchar(inFileName[i]) - 4);
	dir.create(file.path(outDir, outStrDir), showWarnings = FALSE);

	#_______________________________________________________________________________
	## rename chains of domains
	## rename chain of N-terminal domain (first residue set to 1, min is 1)
	pdb$atom$chain[pdb$atom$resno %in% c(1:47)] = "N";
	## create PDB of N-terminal domain
	chainN = atom.select(pdb, resno = c(1:47));
    pdbChainN = trim.pdb(pdb, chainN);

	## rename chain of A domain
	pdb$atom$chain[pdb$atom$resno %in% c(48:116, 221:388)] = "A";
	## create PDB of A domain
	chainA = atom.select(pdb, resno = c(48:116, 221:388));
    pdbChainA = trim.pdb(pdb, chainA);

	## rename chain of B domain
	pdb$atom$chain[pdb$atom$resno %in% c(117:220)] = "B";
	## create PDB of B domain
	chainB = atom.select(pdb, resno = c(117:220));
    pdbChainB = trim.pdb(pdb, chainB);

	## rename chain of C domain (last residue set to 550, max is 533)
	pdb$atom$chain[pdb$atom$resno %in% c(389:550)] = "C";
	## create PDB of C domain
	chainC = atom.select(pdb, resno = c(389:550));
    pdbChainC = trim.pdb(pdb, chainC);

	#_______________________________________________________________________________
	## write out complete PDB structure: all domains with renamed chain(s)
	outFileNameComplete = paste("pdb", ".complete", ".pdb", sep = "");
	write.pdb(pdb, file = paste(outDir, outStrDir, outFileNameComplete, sep = "/"));

	#_______________________________________________________________________________
	## write out single domains
	## N-terminal domain
	outFileNameDom = list();
	doms = list(pdbChainN, pdbChainA, pdbChainB, pdbChainC);
	domsName = c("Ndom", "Adom", "Bdom", "Cdom");

	outFileNameDom[[1]] = paste(domsName[1], "pdb", sep = ".");
	write.pdb(doms[[1]], file = paste(outDir, outStrDir, outFileNameDom[1], sep = "/"));

	## A domain
	outFileNameDom[[2]] = paste(domsName[2], "pdb", sep = ".");
	write.pdb(doms[[2]], file = paste(outDir, outStrDir, outFileNameDom[[2]], sep = "/"));

	## B domain
	outFileNameDom[[3]] = paste(domsName[3], "pdb", sep = ".");
	write.pdb(doms[[3]], file = paste(outDir, outStrDir, outFileNameDom[[3]], sep = "/"));

	## C domain
	outFileNameDom[[4]] = paste(domsName[4], "pdb", sep = ".");
	write.pdb(doms[[4]], file = paste(outDir, outStrDir, outFileNameDom[[4]], sep = "/"));

	#_______________________________________________________________________________
	## write out all domain pairs
	pair.cbn = combn(1:4, 2, simplify = FALSE);
	outFileNameDomPair = list();

	for (j in 1:length(pair.cbn)) {
		## concatenate domains
		domPair = cat.pdb(doms[[pair.cbn[[j]][1]]], doms[[pair.cbn[[j]][2]]],
							renumber = TRUE, rechain = FALSE);
		outFileNameDomPair[[j]] = paste(
						domsName[pair.cbn[[j]][1]],
						domsName[pair.cbn[[j]][2]],
						"pdb", sep = ".");
		write.pdb(domPair, file = paste(outDir, outStrDir, outFileNameDomPair[[j]], sep = "/"));
	}

	#_______________________________________________________________________________
	## POPS all of the above structures (complete, single domain, paired domains)
	## all filenames
	outFileNameAll = list(outFileNameComplete, outFileNameDom, outFileNameDomPair);
	## simplify list structure
	outFileNameAll.l = as.list(unlist(outFileNameAll));
	outFileNameDirAll.l = as.list(paste(outDir, outStrDir, unlist(outFileNameAll), sep = "/"));

	#_______________________________________________________________________________
	## parallelisation
	## initiate cluster for parallel computation 
	#clu = makeCluster(nCore);
	## make parallel functions see predefined variables
	#clusterExport(clu, c("outFileNameDirAll.l"));

	#_______________________________________________________________________________
	## POPS all derived structures of one input structure
	#parLapply(clu, outFileNameDirAll.l, function(x) {
	#	print(paste("./pops --pdb ", x, " --popsOut ", x, ".out", " --popsbOut ", x, ".bout", " --sigmaOut ", x, ".sout", sep = "")); });
	lapply(outFileNameDirAll.l, function(x) {
		system(paste("./pops --pdb ", x, " --popsOut ", x, ".out", " --popsbOut ", x, ".bout", " --sigmaOut ", x, ".sout", " --coarse --rout", sep = "")); });
		#print(paste("./pops --pdb ", x, " --popsOut ", x, ".out", " --popsbOut ", x, ".bout", " --sigmaOut ", x, ".sout", " --coarse --rout", sep = "")); });

	## release memory of parallelised structure
	#stopCluster(clu);

	#_______________________________________________________________________________
	## POPSCOMP : compute $DeltaSASA in complex [(X + Y) - XY]
	## for all domain pairs
	outFileNameDeltaSasa = list();
	for (j in 1:length(pair.cbn)) {
		domXY.molSasa.outFileName = paste(outFileNameDomPair[[j]], "out.rpopsMolecule", sep = ".");
		domXY.molSasa = read.table(paste(outDir, outStrDir, domXY.molSasa.outFileName, sep = "/"), header = TRUE);
		domX.molSasa.outFileName = paste(outFileNameDom[[pair.cbn[[j]][1]]], "out.rpopsMolecule", sep = ".");
		domX.molSasa = read.table(paste(outDir, outStrDir, domX.molSasa.outFileName, sep = "/"), header = TRUE);
		domY.molSasa.outFileName = paste(outFileNameDom[[pair.cbn[[j]][2]]], "out.rpopsMolecule", sep = ".");
		domY.molSasa = read.table(paste(outDir, outStrDir, domY.molSasa.outFileName, sep = "/"), header = TRUE);
		## deltaSASA
		delta.molSasa = (domX.molSasa + domY.molSasa) - domXY.molSasa; 
		outFileNameDeltaSasa[[j]] = paste(
						domsName[pair.cbn[[j]][1]],
						domsName[pair.cbn[[j]][2]],
						"dsasaMolecule", sep = ".");
		write.table(delta.molSasa, file = paste(outDir, outStrDir, outFileNameDeltaSasa[[j]], sep = "/"));
	}
}

#===============================================================================
