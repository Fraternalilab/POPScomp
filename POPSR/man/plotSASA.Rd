\name{plotSASAatom}
\alias{plotSASAatom}
\alias{plotSASAresidue}
\alias{plotSASAchain}
\title{Plot SASA values of a POPS run}
\usage{
plotSASAatom(sasaFile = "id.rpopsAtom", plotFile = "sasa_atom.png")

plotSASAresidue(sasaFile = "id.rpopsResidue", plotFile = "sasa_residue.png")

plotSASAchain(sasaFile = "id.rpopsChain", plotFile = "sasa_chain.png")
}
\arguments{
\item{sasaFile}{R-format SASA table written by POPS with the \code{--rout}
  switch; its name is \code{<routPrefix>.rpopsAtom}, \code{.rpopsResidue} or
  \code{.rpopsChain}, where the default prefix of POPS is \code{id}.}

\item{plotFile}{name of the PNG file to create.}
}
\value{
The name of the written plot file, invisibly; the plot is written to
\code{plotFile}.
}
\description{
Plot the solvent accessible surface area per atom, per residue (total,
hydrophobic and hydrophilic) or per chain.
}
\examples{
\dontrun{
plotSASAresidue("id.rpopsResidue", "sasa_residue.png")
}
}
