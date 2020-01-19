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
          p("The POPS program computes the Solvent Accessible Surface Area (SASA)
            of a given PDB structure. If the structure was composed of more than one chain
            containing protein or RNA/DNA, POPScomp creates internally all pair combinations
            of chains to compute the buried SASA upon complexation. Details of the procedure
            are explained in the published papers on POPS and POPSCOMP."),
          br(),
          p("POPScomp shows tabs for atom, residue, chain and molecule SASAs.
            The tables are initialised without any values and therefore the user sees
            the table header and below the notice 'Showing 0 to 0 of 0 entries'.
            After selecting a PDB file and pressing 'run POPScomp', the sever runs
            the POPS program on components of the PDB file. Because that computation
            is a system call, the success of the computation is returned
            as exit code and shown below the 'run POPScomp' button:"
          ),
          p("* 0 - Success"),
          p("* 1 - Catchall for general errors"),
          p("* 2 - Misuse of shell builtins (according to Bash documentation"),
          p("* 126 - Command invoked cannot execute"),
          p("* 127 - Command not found"),
          p("* 128 - Invalid argument to exit"),
          p("* 128+n - Fatal error signal 'n'"),
          p("* 130 - Script terminated by Control-C"),
          p("* 255* - Exit status out of range")
        ),
        tabPanel("About",
            h5("This is version 3.0.4 of the", a("POPScomp server",
                                    href="http://popscom.org:3838")),
            p("The POPScomp server is based on two software packages:"),
            p("1. A GNU Autotools package of the POPS C program."),
            p("2. An R package containing a Shiny server",
              "to interface the POPS program and POPSCOMP functionality."),
            p("Since April 2019, the POPS program (POPSC) and the",
              "POPScomp Shiny server (POPSR) are being co-developed."),
            h5("Source code and detailed information can be found on
               Fraternali Lab's", a("POPScomp GitHub page",
                                    href="https://github.com/Fraternalilab/POPScomp")),
            br(),
            h5("POPScomp is part of the FunPDBe", a("FunPDBe resources",
                                    href="https://www.ebi.ac.uk/pdbe/funpdbe/deposition")),
            br(),
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
