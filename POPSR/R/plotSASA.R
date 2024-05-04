#! /usr/bin/R
#===============================================================================
# POPSR package
# Plot SASA
# (C) 2019-2023 Jens Kleinjung and Franca Fraternali
#===============================================================================

#_______________________________________________________________________________
# plot SASA per atom
plotSASAatom = function() {
	sasa.atom = read.table("rpopsAtom_pops.out", header = TRUE);

	png("sasa_atom.png");

	g = ggplot(sasa.atom, aes(x = AtomNr, y = SASA.A.2)) + 
			geom_point();
	g = g + xlab("atom number") +
			ylab(expression(paste("SASA / ", A^2, sep = "")));
	g = g + scale_alpha(guide = 'none');
	g = g + theme(plot.margin = unit(c(1,1,1,1), "cm"));
	g = g + theme(axis.text = element_text(size = 14),
			axis.title = element_text(size = 14));
	plot(g);

	dev.off();
}

#_______________________________________________________________________________
# plot SASA per residue 
plotSASAresidue = function() {
	sasa.residue = read.table("rpopsResidue_pops.out", header = TRUE);

	png("sasa_residue.png");

	g = ggplot(data = sasa.residue, aes(x = ResidNr, y = Total.A.2)) + 
			geom_point();
	g = g + geom_point(data = sasa.residue, aes(x = ResidNr, y = Phob.A.2), colour = "green");
	g = g + geom_point(data = sasa.residue, aes(x = ResidNr, y = Phil.A.2), colour = "blue");
	g = g + xlab("residue number") +
		ylab(expression(paste("SASA / ", A^2, sep = "")));
	g = g + theme(plot.margin = unit(c(1,1,1,1), "cm"));
	g = g + theme(axis.text = element_text(size = 14),
			axis.title = element_text(size = 14));

	plot(g);

	dev.off();
}

#_______________________________________________________________________________
# plot SASA per chain and molecule
plotSASAchain = function() {
	sasa.chain = read.table("rpopsChain_pops.out", header = TRUE);

	colours = rep(c("green", "blue", "black"), dim(sasa.chain)[1]);

	chain.id = rep(sasa.chain$Id, 3);
	chain.sasa = c(sasa.chain$Phob.A.2, sasa.chain$Phil.A.2, sasa.chain$Total.A.2);
	chain.sasatype = factor(rep(c("phob", "phil", "total"), each = 2));
	chain.df = data.frame(chain.id, chain.sasa, chain.sasatype, paste(chain.id, chain.sasatype, sep = "_"));
	colnames(chain.df) = c("id", "sasa", "sasatype", "id.sasatype");

	png("sasa_chain.png");

	g = ggplot(data = chain.df, aes(x = id.sasatype, y = sasa, fill = sasatype)) + 
			geom_col(colour = "black");
	g = g + scale_fill_manual("", values = colours, drop = FALSE);

	g = g + xlab("chain") +
			ylab(expression(paste("SASA / ", A^2, sep = "")));
	g = g + theme(plot.margin = unit(c(1,1,1,1), "cm"));
	g = g + theme(axis.text = element_text(size = 14),
			axis.title = element_text(size = 14));

	plot(g);

	dev.off();
}

#===============================================================================
