#! /usr/bin/R
#===============================================================================
# plot SASA per residue 
#===============================================================================

library("ggplot2");

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

#===============================================================================
