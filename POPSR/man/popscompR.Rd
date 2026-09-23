\name{popscompR}
\alias{popscompR}
\alias{POPSR}
\alias{POPScomp}
\title{POPScomp: buried surface areas of the chains and chain pairs of a complex}
\usage{
popscompR(inputPDB, outDir, coarse = FALSE)
}
\arguments{
\item{inputPDB}{name of the input PDB file, relative to \code{outDir}.}

\item{outDir}{directory holding the input structure; all output files of the
  run are written there.}

\item{coarse}{if \code{TRUE}, POPS computes one centre per residue
  (\code{--coarse}) instead of one per atom; the atom level is then absent
  from the results.}
}
\value{
A list of the difference (buried) SASA tables, one element per resolution
level (\code{rpopsAtom}, \code{rpopsResidue}, \code{rpopsChain},
\code{rpopsMolecule}; without \code{rpopsAtom} in coarse mode), each holding
one table per chain pair. The same tables are written to \code{outDir} as
\code{deltaSASA.<level>}. Returns \code{0} for structures with fewer than two
usable chains, because POPScomp processes complexes.
}
\description{
The input structure is split into its chains. POPS is run on every isolated
chain and on every pair of chains, and the buried surface area of each chain
pair is the SASA of the two isolated chains minus the SASA of the pair. The
POPS program is located through the \code{POPS_BIN} environment variable, the
\code{PATH}, or the installation path of the Shiny image, in that order.
}
\examples{
\dontrun{
## 'outDir' holds the input structure and receives all output files
popscompR("1f3r.pdb", outDir = tempdir())
}
}
