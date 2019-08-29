#! /usr/bin/R
#===============================================================================
## Shiny application for R interface of POPS
# (C) 2019 Jens Kleinjung and Franca Fraternali
#===============================================================================

library("shiny")
library("DT")
library("digest")
library("RPOPS")
library("bio3d")
library("parallel")

#_______________________________________________________________________________
#_______________________________________________________________________________
# POPS UI
ui <- fluidPage(

  #_______________________________________________________________________________
  titlePanel("POPScomp", windowTitle = "POPScomp"),

  #_______________________________________________________________________________
  ## sidebar layout
  sidebarLayout(
    sidebarPanel(
      #_______________________________________________________________________________
      ## input 1 : PDB file
      ## PDB identifier
      textInput(inputId = "pdbentry",
                label = "Enter PDB entry:",
                value = ""),
      ## PDB upload
      fileInput("PDBfile", "OR upload PDB file",
                multiple = FALSE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".pdb")),

      ## horizontal line
      tags$hr(),

      #_______________________________________________________________________________
      ## input 2 : POPS mode
      selectInput(inputId = "popsmode",
                label = "Resolution:",
                choices = c("atomistic", "coarse")),
      ## input 3 : Probe (= solvent) radius
      numericInput(inputId = "rprobe",
                   label = "Solvent radius [Angstrom]:",
                   value = 1.4),

      ## horizontal line
      tags$hr(),

      #_______________________________________________________________________________
      ## input 4 : Action button
      actionButton("popscomp", label = "run POPScomp", class = "btn-primary"),

      textOutput("nil"),

      ## horizontal line
      tags$hr()
    ),

    #_______________________________________________________________________________
    ## main panel for output
    mainPanel(
      tabsetPanel(
        #_______________________________________________________________________________
        ## Atom SASA table
        tabPanel("Atom",
          DT::dataTableOutput("popsAtom")
        ),
        #_______________________________________________________________________________
        ## Residue SASA table
        tabPanel("Residue",
          DT::dataTableOutput("popsResidue")
        ),
        #_______________________________________________________________________________
        ## Chain SASA table
        tabPanel("Chain",
          DT::dataTableOutput("popsChain")
        ),
        #_______________________________________________________________________________
        ## Molecule SASA table
        tabPanel("Molecule",
          DT::dataTableOutput("popsMolecule")
        ),
        #_______________________________________________________________________________
        ## Summary SASA table
        tabPanel("Summaries",
                 p("Summaries of SASAs")
        ),
        #_______________________________________________________________________________
        ## Plots
        tabPanel("Plots",
          p("Plots of SASAs")
        ),
        #_______________________________________________________________________________
        ## Readme
        tabPanel("Readme",
          p("The POPS program computes the Solvent Accessible Surface Area (SASA)
            of a given PDB structure. If the structure is composed of more than one chain
            containing protein or RNA/DNA, POPScomp creates internally all pair combinations
            of chains to compute the buried SASA upon complexation. Details of those routines
            are explained in the published papers on POPS and POPSCOMP."),
          br(),
          p("POPScomp shows tabs for atom, residue, chain and molecule SASAs.
            The tables are initialised without any values and therefore the user sees
            the table header and below the notice 'No data available in table'.
            After selecting a PDB identifier or file and pressing 'run POPScomp',
			the sever runs the POPS program on components of the PDB file
			and the tables automatically refresh to show the resulting SASA values.
			Because running POPS is a system call, the success of the computation
			is returned as exit code and shown below the 'run POPScomp' button:"
          ),
          p("* 0 - Success"),
          p("* 1 - Catchall for general errors"),
          p("* 2 - Misuse of shell builtins (according to Bash documentation)"),
          p("* 126 - Command invoked cannot execute"),
          p("* 127 - Command not found"),
          p("* 128 - Invalid argument to exit"),
          p("* 128+n - Fatal error signal 'n'"),
          p("* 130 - Script terminated by Control-C"),
          p("* 255* - Exit status out of range"),
		  br(),
		  p("Results will be stored on the server for maximally one day.
			Upon high demand the storage time might be reduced 30 minutes.
			Please download your results via the 'Download' buttons.")
        ),
		  #_______________________________________________________________________________
		  ## About
		  tabPanel("About",
          h5("This is version 3.0.0 of the", a("POPScomp server",
              href="http://popscom.org:3838")),
          p("The server automatically recognises PDB identifiers and multi-chain structures.
              Output comprises downloadable SASA tables and graphs shown on the Shiny server pages."),
          br(),
          p("The POPScomp server is based on two software packages:"),
          p("1. A GNU Autotools package of the POPS C program that computes SASA for a given structure."),
          p("2. An R package containing an R program that a) splits complexes into single and pair
              components to compute buried SASA using POPSC and b) provides a Shiny server to interface the R program."),
          p("Since April 2019 the POPS program (POPSC) and the",
              "POPScomp Shiny server (POPSR) are being co-developed."),
          br(),
          h5("Source code and detailed information can be found on
              Fraternali Lab's", a("POPScomp GitHub page",
              href="https://github.com/Fraternalilab/POPScomp")),
          p("Please use that site for bug reports and other comments."),
          br(),
          h5("POPScomp is part of the", a("FunPDBe resources",
              href="https://www.ebi.ac.uk/pdbe/funpdbe/deposition")),
          br(),
          p("Usage of the server is free, the code license is GPL3."),
          p("Authors:",
              "Franca Fraternali (franca.fraternali@kcl.ac.uk)",
              "and Jens Kleinjung (jens@jkleinj.eu)")
        )
      )
    )
  )
)

#_______________________________________________________________________________
# server routines
server <- function(input, output) {

  #_______________________________________________________________________________
  ## output 1 : display input PDB entry
  output$pdbentry <- renderText({
    input$pdbentry
  })

  #_______________________________________________________________________________
  ## output 2 : display input POPS mode
  output$popsmode <- renderText({
    input$popsmode
  })

  #_______________________________________________________________________________
  ## output 3 : display input probe radius
  output$rprobe <- renderText({
    input$rprobe
  })

  #_______________________________________________________________________________
  ## silent reactive output 4: download PDB entry or upload input file
  ## run popscompR on specified PDB file
  ## Comment: The 'fileInput' function returns the object 'input$file1',
  ##  a list of four elements, of which the fourth element is the
  ##  path to the temporary file. See below 'input$PDBfile[[4]]'.
  output$nil <- eventReactive(input$popscomp, {
    ## to proceed, we require one PDB identifier or uploaded PDB file
    validate(need(((input$pdbentry != "") || (! is.null(input$PDBfile))),
          message = "No PDB source input!"))
    ## to proceed, we need one unspecified PDB identifier
    validate(need(((input$pdbentry == "") || (is.null(input$PDBfile))),
          message = "Two PDB sources input!"))

    ## set/create output directory
    outdir = tempdir();
    setwd(outdir);

    ## download/copy PDB structure
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

    #_______________________________________________________________________________
    ## run popscompR
    sasa.l = popscompR(inputPDB, outdir);

  })

  #_______________________________________________________________________________
  ## output 5 : display POPS output for atoms, residues, chains and molecule
  #_____________________________________
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
  atomOutput = reactiveValues(highlight = NULL, data = NULL)
  atomOutputData = sasa.l[[1]][[1]];
  ## render output data as table
  output$popsAtom = DT::renderDataTable({
    atomOutput$data = atomOutputData
  })

  #_____________________________________
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
  residueOutput = reactiveValues(highlight = NULL, data = NULL)
  residueOutputData = sasa.l[[1]][[1]];
  ## render output data as table
  output$popsResidue = DT::renderDataTable({
    residueOutput$data = residueOutputData
  })

  #_____________________________________
  ## chain
  chain_null.df = data.frame(
                    ChainNr = integer(),
                    ChainNe = character()
  )
  chainOutput = reactiveValues(highlight = NULL, data = NULL)
  chainOutputData = sasa.l[[1]][[1]];
  ## render output data as table
  output$popsChain = DT::renderDataTable({
    chainOutput$data = chainOutputData
  })

  #_____________________________________
  ## molecule
  molecule_null.df = data.frame(
                      Phob = double(),
                      Phil = double(),
                      Total = double()
  )
  moleculeOutput = reactiveValues(highlight = NULL, data = NULL)
  moleculeOutputData = sasa.l[[1]][[1]];
  ## render output data as table
  output$popsMolecule = DT::renderDataTable({
    moleculeOutput$data = moleculeOutputData
  })
}

#_______________________________________________________________________________
# run the Shiny app
shinyApp(ui, server)

#===============================================================================
