#! /usr/bin/R
#===============================================================================
# plot SASA per chain and molecule
#===============================================================================

library("ggplot2");

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

#===============================================================================
