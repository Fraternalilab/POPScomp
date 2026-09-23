\name{runPOPSR}
\alias{runPOPSR}
\title{Run the POPScomp Shiny application}
\usage{
runPOPSR()
}
\value{
Called for its side effect: the Shiny application is started and does not
return until it is closed.
}
\description{
Start the Shiny interface of POPScomp, which runs POPS on an uploaded or
downloaded structure and shows the SASA and buried-SASA tables. The
application additionally needs the packages \code{shiny}, \code{DT},
\code{digest}, \code{markdown} and \code{shinysky}.
}
\examples{
\dontrun{
runPOPSR()
}
}
