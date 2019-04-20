#! /usr/bin/R
#===============================================================================
## shiny application for R interface of POPS
# (C) 2019 Jens Kleinjung
#===============================================================================

library(shiny)
library(bio3d)
library(DT)
library(digest)

#_______________________________________________________________________________
# POPS UI
ui <- fluidPage(

  titlePanel("POPScomp (version 3.0.0)", windowTitle = "POPScomp"),

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
  ## Comment: The 'fileInput' function returns the object 'input$file1',
  ##  a list of four elements, of which the fourth element is the
  ##  path to the temporary file.
  output$nil <- eventReactive(input$popscomp, {
    ## to proceed, we require one PDB identifier or uploaded PDB file
    validate(need(((input$pdbentry != "") || (! is.null(input$PDBfile))),
          message = "No PDB source input!"))
    ## to proceed, we need one unspecified PDB identifier
    validate(need(((input$pdbentry == "") || (is.null(input$PDBfile))),
          message = "Two PDB sources input!"))

    ## set output directory
    mainDir = "/tmp"
    ## creates random string based on subsecond time
    subDir = paste0("POPScomp_", digest(format(Sys.time(), "%H:%M:%OS3")))
    outdir = paste(mainDir, subDir, sep = "/")
    dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
    setwd(file.path(mainDir, subDir))
    system("ln -s ~/0/POPS/src/pops")

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

    ## run POPS as system command
    command = paste("./pops --outDirName", outdir,
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
