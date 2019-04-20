#! /usr/bin/R
#===============================================================================
# plot SASA per atom
#===============================================================================

library("ggplot2");

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

#===============================================================================
