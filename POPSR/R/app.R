#! /usr/bin/R
#===============================================================================
# Shiny application as interface of POPS & POPSCOMP
# (C) 2019 -2020 Jens Kleinjung and Franca Fraternali
#===============================================================================

library(shiny)
library(bio3d)
library(DT)
library(digest)

#_______________________________________________________________________________
# POPS UI
ui <- fluidPage(

  titlePanel("POPScomp", windowTitle = "POPScomp"),

  ## sidebar layout
  sidebarLayout(
    sidebarPanel(
      ## i1.1
      textInput(inputId = "pdbentry",
                label = "Enter PDB entry:",
                value = ""),
      ## i1.2
      fileInput("PDBfile", "OR upload PDB file",
                multiple = FALSE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".pdb")),

      # horizontal line
      tags$hr(),

      ## i2
      selectInput(inputId = "popsmode",
                label = "Resolution:",
                choices = c("atomistic", "coarse")),
      ## i3
      numericInput(inputId = "rprobe",
                   label = "Solvent radius [Angstrom]:",
                   value = 1.4),

      # horizontal line
      tags$hr(),

      # i4 action button
      actionButton("popscomp", label = "run POPScomp", class = "btn-primary"),

      textOutput("nil"),

      # horizontal line
      tags$hr()
    ),

    ## main panel for output
    mainPanel(
      tabsetPanel(
        tabPanel("Atom",
          DT::dataTableOutput("popsAtom")
        ),
        tabPanel("Residue",
          DT::dataTableOutput("popsResidue")
        ),
        tabPanel("Chain",
          DT::dataTableOutput("popsChain")
        ),
        tabPanel("Molecule",
          DT::dataTableOutput("popsMolecule")
        ),
        tabPanel("Readme",
		      h3("Method"),
          p("The POPScomp server invokes the POPS program to compute the
		      Solvent Accessible Surface Area (SASA) of a given PDB structure.
			    This is the POPS functionality, which has been currently implemented.
			    Soon to come: The POPSCOMP functionality for protein or RNA/DNA complexes,
			    where the POPScomp server creates internally
			    all pair combinations of chains to compute the buried SASA upon complexation.
			    Details of those functionalities are explained in the published papers
			    on POPS and POPSCOMP, see 'About' tab for the list of publications."),
		      h3("Results"),
          p("The POPScomp Shiny app shows tabs for atom, residue, chain and molecule SASAs.
            The tables are initialised without any values; therefore, before 'run POPScomp' execution,
            the user sees only the table header and below the notice 'Showing 0 to 0 of 0 entries'.
            After selecting a PDB file and pressing 'run POPScomp', the server runs
            the POPS program on the selected PDB file. The output are SASA tables, which are
            automatically loaded into the respective tabs.
            Because that computation is a shell command, the success of the computation is
            returned as exit code and shown below the 'run POPScomp' button.
            See 'Exit Codes' tab for details."
          ),
		      h3("Help"),
          p("In case the server does not work as expected or server-related issues
		        need clarification, please email the maintainers:
			      Jens Kleinjung (jens@jkleinj.eu) and
            Franca Fraternali (franca.fraternali@kcl.ac.uk).")
        ),
        tabPanel("About",
			    h3("Server"),
			    h4("Software"),
          h5("This is version 3.0.4 of the", a("POPScomp server",
                                href="http://popscom.org:3838")),
          p("The POPScomp server is based on two software packages:"),
          p("1. A GNU Autotools package of the POPS C program."),
          p("2. An R package containing this Shiny server
             to interface the POPS program and POPSCOMP functionality."),
          p("Since April 2019, the POPS program (POPSC) and the
             POPScomp Shiny server (POPSR) are being co-developed."),
			    h5("Source code and detailed information can be found on
              Fraternali Lab's ", a("POPScomp GitHub page",
                                href="https://github.com/Fraternalilab/POPScomp")),
			    h4("EBI PDBe-KB"),
          h5("POPScomp is part of the ", a("FunPDBe resources",
                                href="https://www.ebi.ac.uk/pdbe/funpdbe/deposition")),
			    h3("References"),
			    h5("Users publishing results obtained with the program and
			        its applications should acknowledge its use by citation."),
			    h4("Implicit solvent"),
			    h5("Fraternali, F. and van Gunsteren, W.F.
			        An efficient mean solvation force model for use in
			        molecular dynamics simulations of proteins in aqueous solution.
			        Journal of Molecular Biology 256 (1996) 939-948."),
			    h4("POPS method"),
			    h5("Fraternali, F. and Cavallo, L.
			        Parameter optimized surfaces (POPS): analysis of key interactions
			        and conformational changes in the ribosome.
			        Nucleic Acids Research 30 (2002) 2950-2960."),
			    h4("POPS server"),
			    h5("Cavallo, L., Kleinjung, J. and Fraternali, F.
			        POPS: A fast algorithm for solvent accessible surface areas
			        at atomic and residue level.
			        Nucleic Acids Research 31 (2003) 3364-3366."),
			    h4("POPSCOMP server"),
			    h5("Kleinjung, J. and Fraternali, F.
			        POPSCOMP: an automated interaction analysis of biomolecular complexes.
			        Nucleic Acids Research 33 (2005) W342-W346."),
			    h3("License"),
			    h5("Usage of the software and server is free under the
			        GNU General Public License v3.0."),
			    h4("Copyright Holders, Authors and Maintainers"),
			    h5("2002-2020 Franca Fraternali (author, maintainer)"),
			    h5("2008-2020 Jens Kleinjung (author, maintainer)"),
			    h4("Contributors"),
			    h5("2002 Kuang Lin and Valerie Hindie (translation to C)"),
			    h5("2002 Luigi Cavallo (parametrisation)")
        ),
        tabPanel("Exit Codes",
          h3("Shiny exit codes"),
          p("* No PDB source input! - Enter PDB identifier or upload PDB file from local file system
            at the top of the side panel."),
          p("* Two PDB sources input! - Only one PDB source is accepted per computation. Refresh the
            browser page and either speficy a PDB identifier or upload a PDB file, not both."),
          h3("Shell command exit codes"),
          p("* 0 - Success"),
          p("* 1 - Catchall for general errors"),
          p("* 2 - Misuse of shell builtins (according to Bash documentation"),
          p("* 126 - Command invoked cannot execute"),
          p("* 127 - Command not found"),
          p("* 128 - Invalid argument to exit"),
          p("* 128+n - Fatal error signal 'n'"),
          p("* 130 - Script terminated by Control-C"),
          p("* 255* - Exit status out of range")
        )
      )
    )
  )
)

#_______________________________________________________________________________
# server routines
server <- function(input, output) {

  ## o1.1 display input PDB entry
  output$pdbentry <- renderText({
    input$pdbentry
  })

  ## o2 display input POPS mode
  output$popsmode <- renderText({
    input$popsmode
  })

  ## o3 display input probe radius
  output$rprobe <- renderText({
    input$rprobe
  })

  ## o4 download PDB entry or upload input file
  ## run POPS on specified PDB file
  ## Comments:
  ## - Expects 'pops' binary as /usr/local/bin/pops.
  ## - The App will set its own working directory to a temporary directory,
  ## to which the specified PDB file will be up/down-loaded.
  ## - For uploaded files: The 'fileInput' function returns the object 'input$file1',
  ##  a list of four elements, of which the fourth element contains the
  ##  path to the temporary file.
  output$nil <- eventReactive(input$popscomp, {
    ## to proceed, we require one PDB identifier or uploaded PDB file
    validate(need(((input$pdbentry != "") || (! is.null(input$PDBfile))),
          message = "No PDB source input!"))
    ## to proceed, we need one unspecified PDB identifier
    validate(need(((input$pdbentry == "") || (is.null(input$PDBfile))),
          message = "Two PDB sources input!"))

    mainDir = "/tmp"
    ## creates random string based on subsecond time
    subDir = paste0("POPScomp_", digest(format(Sys.time(), "%H:%M:%OS3")))
    outdir = paste(mainDir, subDir, sep = "/")
    dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
    setwd(file.path(mainDir, subDir))

    ## download (PDB database) or upload (local file system) the PDB structure
    if (input$pdbentry != "") {
      ## get.pdb downloads the PDB structure from the database
      get.pdb(input$pdbentry, path = outdir);
      inputPDB = paste(input$pdbentry, "pdb", sep = ".")
    } else {
      ## move uploaded PDB file from its temporary directory to output directory
      system(paste("mv ", input$PDBfile[[4]], " ", outdir, "/",
                   input$PDBfile[[1]],  sep = ""))
      inputPDB = input$PDBfile[[1]]
    }

    ## run POPS as system command
    command = paste("/usr/local/bin/pops --outDirName", outdir,
                    "--rout --atomOut --residueOut --chainOut",
                    "--pdb", inputPDB, "1> POPScomp.o 2> POPScomp.e");
    system_status = system(command)
    paste("Exit code:", system_status)
  })

  ## o5 display POPS output for atoms, residues, chains and molecule
  ## atom
  ## empty dataframe with column names
  ## that wil show up as empty table before POPS has been run
  atom_null.df = data.frame(
                      AtomNr = integer(),
                      AtomNe = character(),
                      Residue = character(),
                      Chain = character(),
                      ResidNr = integer(),
                      iCode = integer(),
                      SASA = double(),
                      Q = double(),
                      N = integer(),
                      AtomTp = integer(),
                      AtomGp = integer(),
                      Surf = double()
  )
  write.table(atom_null.df, file = "pops.out.rpopsAtom")
  ## reactive data: update output when file content changes
  atomOutput = reactiveValues(highlight = NULL, data = NULL)
  atomOutputData = reactiveFileReader(2000, NULL, "pops.out.rpopsAtom",
                                         read.table, header = TRUE)
  ## render output data as table
  output$popsAtom = DT::renderDataTable({
    atomOutput$data = atomOutputData()
  })

  ## residue
  ## create empty dataframe in r with column names
  residue_null.df = data.frame(
                      ResidNe = character(),
                      Chain = character(),
                      ResidNr = integer(),
                      iCode = integer(),
                      Phob = double(),
                      Phil = double(),
                      Total = double(),
                      Q = double(),
                      N = integer(),
                      Surf = double()
  )
  write.table(residue_null.df, file = "pops.out.rpopsResidue")
  residueOutput = reactiveValues(highlight = NULL, data = NULL)
  residueOutputData = reactiveFileReader(2000, NULL, "pops.out.rpopsResidue",
                                      read.table, header = TRUE)
  output$popsResidue = DT::renderDataTable({
    residueOutput$data = residueOutputData()
  })

  ## chain
  chain_null.df = data.frame(
                    ChainNr = integer(),
                    ChainNe = character()
  )
  write.table(chain_null.df, file = "pops.out.rpopsChain")
  chainOutput = reactiveValues(highlight = NULL, data = NULL)
  chainOutputData = reactiveFileReader(2000, NULL, "pops.out.rpopsChain",
                                      read.table, header = TRUE)
  output$popsChain = DT::renderDataTable({
    chainOutput$data = chainOutputData()
  })

  ## molecule
  molecule_null.df = data.frame(
                      Phob = double(),
                      Phil = double(),
                      Total = double()
  )
  write.table(molecule_null.df, file = "pops.out.rpopsMolecule")
  moleculeOutput = reactiveValues(highlight = NULL, data = NULL)
  moleculeOutputData = reactiveFileReader(2000, NULL, "pops.out.rpopsMolecule",
                                      read.table, header = TRUE)
  output$popsMolecule = DT::renderDataTable({
    moleculeOutput$data = moleculeOutputData()
  })
}

#_______________________________________________________________________________
# run the Shiny app
shinyApp(ui, server)

#===============================================================================
