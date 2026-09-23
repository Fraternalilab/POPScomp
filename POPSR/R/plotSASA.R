#! /usr/bin/R
#===============================================================================
# POPSR package
# Plot SASA
# The input files are the R-format tables written by POPS with the '--rout'
#   switch, whose names are '<routPrefix>.rpopsAtom', '.rpopsResidue' and
#   '.rpopsChain'; the default prefix of POPS is 'id'.
# (C) 2019-2026 Jens Kleinjung and Franca Fraternali
#===============================================================================

## column names of the POPS tables, used in the ggplot2 aesthetics below
utils::globalVariables(c("AtomNr", "SASA.A.2", "ResidNr", "sasa", "type",
                         "id", "sasatype"))

#_______________________________________________________________________________
# plot SASA per atom
plotSASAatom = function(sasaFile = "id.rpopsAtom", plotFile = "sasa_atom.png") {
	sasa.atom = read.table(sasaFile, header = TRUE);

	png(plotFile);
	on.exit(dev.off());

	g = ggplot2::ggplot(sasa.atom, ggplot2::aes(x = AtomNr, y = SASA.A.2)) + 
			ggplot2::geom_point();
	g = g + ggplot2::xlab("atom number") +
			ggplot2::ylab(expression(paste("SASA / ", ring(A)^2, sep = "")));
	g = g + ggplot2::scale_alpha(guide = 'none');
	g = g + ggplot2::theme(plot.margin = ggplot2::unit(c(1,1,1,1), "cm"));
	g = g + ggplot2::theme(axis.text = ggplot2::element_text(size = 14),
			axis.title = ggplot2::element_text(size = 14));
	plot(g);
}

#_______________________________________________________________________________
# plot SASA per residue 
plotSASAresidue = function(sasaFile = "id.rpopsResidue", plotFile = "sasa_residue.png") {
	sasa.residue = read.table(sasaFile, header = TRUE);

	png(plotFile);
	on.exit(dev.off());

	## total, hydrophobic and hydrophilic SASA of each residue;
	## the colours are distinguishable with colour vision deficiency
	sasa.long = data.frame(
		ResidNr = rep(sasa.residue$ResidNr, 3),
		sasa = c(sasa.residue$SASA.A.2, sasa.residue$Phob.A.2, sasa.residue$Phil.A.2),
		type = factor(rep(c("total", "hydrophobic", "hydrophilic"),
				each = nrow(sasa.residue)),
			levels = c("total", "hydrophobic", "hydrophilic")));

	g = ggplot2::ggplot(sasa.long, ggplot2::aes(x = ResidNr, y = sasa, colour = type)) +
			ggplot2::geom_point();
	g = g + ggplot2::scale_colour_manual("",
			values = c(total = "black", hydrophobic = "#E69F00", hydrophilic = "#0072B2"));
	g = g + ggplot2::xlab("residue number") +
		ggplot2::ylab(expression(paste("SASA / ", ring(A)^2, sep = "")));
	g = g + ggplot2::theme(plot.margin = ggplot2::unit(c(1,1,1,1), "cm"));
	g = g + ggplot2::theme(axis.text = ggplot2::element_text(size = 14),
			axis.title = ggplot2::element_text(size = 14));

	plot(g);
}

#_______________________________________________________________________________
# plot SASA per chain and molecule
plotSASAchain = function(sasaFile = "id.rpopsChain", plotFile = "sasa_chain.png") {
	sasa.chain = read.table(sasaFile, header = TRUE, colClasses = c(Id = "character"));

	## one bar per chain and SASA type
	chain.df = data.frame(
		id = rep(as.character(sasa.chain$Id), 3),
		sasa = c(sasa.chain$Phob.A.2, sasa.chain$Phil.A.2, sasa.chain$SASA.A.2),
		sasatype = factor(rep(c("hydrophobic", "hydrophilic", "total"),
				each = nrow(sasa.chain)),
			levels = c("hydrophobic", "hydrophilic", "total")));

	png(plotFile);
	on.exit(dev.off());

	## one group of bars per chain
	g = ggplot2::ggplot(data = chain.df,
			ggplot2::aes(x = id, y = sasa, fill = sasatype)) + 
			ggplot2::geom_col(colour = "black", position = ggplot2::position_dodge());
	## colours distinguishable with colour vision deficiency
	g = g + ggplot2::scale_fill_manual("",
			values = c(hydrophobic = "#E69F00", hydrophilic = "#0072B2", total = "grey70"),
			drop = FALSE);

	g = g + ggplot2::xlab("chain") +
			ggplot2::ylab(expression(paste("SASA / ", ring(A)^2, sep = "")));
	g = g + ggplot2::theme(plot.margin = ggplot2::unit(c(1,1,1,1), "cm"));
	g = g + ggplot2::theme(axis.text = ggplot2::element_text(size = 14),
			axis.title = ggplot2::element_text(size = 14));

	plot(g);
}

#===============================================================================
