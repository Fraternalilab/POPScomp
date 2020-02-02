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

  titlePanel(title=div(img(
                        src="POPScomp.png",
                        width = 150, height = 120,
                        style = "margin 5px 5px"
                      ),
                      "POPScomp")),
  #titlePanel("POPScomp", windowTitle = "POPScomp"),

  ## sidebar layout
  sidebarLayout(
    sidebarPanel(
      ## i1.1
      textInput(inputId = "pdbentry",
                label = "Enter PDB ID (4 characters):",
                value = ""),

      ## i1.2
      fileInput("PDBfile", "OR upload PDB file",
                multiple = FALSE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".pdb")),
      ## horizontal line
      tags$hr(),

      ## i2
      selectInput(inputId = "popsmode",
                label = "Resolution:",
                choices = c("atomistic", "coarse")),

      ## i3
      numericInput(inputId = "rprobe",
                   label = "Solvent radius [Angstrom]:",
                   value = 1.4),
      tags$hr(),

      # i4 action button
      actionButton("popscomp", label = "run POPScomp", class = "btn-primary"),
      textOutput("nil"),
      tags$hr(),

      ## i5 download button
      downloadButton('downloadAllResults', 'Download All Results')
    ),

    ## main panel for output
    mainPanel(
      tabsetPanel(
        tabPanel("Atom",
          tabsetPanel(
            tabPanel("Input Structure",
              DT::dataTableOutput("popsSASAAtom"),
              tags$hr(),
              downloadButton('downloadAtomSASA', 'Download Atom SASA')
            ),
            tabPanel("DeltaSASA",
              DT::dataTableOutput("popsDeltaSASAAtom"),
              tags$hr(),
              downloadButton('downloadAtomDeltaSASA', 'Download Atom DeltaSASA')
            ),
            tabPanel("Isolated Chains",
              DT::dataTableOutput("popsIsoSASAAtom"),
              tags$hr(),
              downloadButton('downloadAtomIsoSASA', 'Download Isolated-Chains Atom SASA')
            )
          )
        ),
        tabPanel("Residue",
          tabsetPanel(
            tabPanel("Input Structure",
              DT::dataTableOutput("popsSASAResidue"),
              tags$hr(),
              downloadButton('downloadResidueSASA', 'Download Residue SASA')
            ),
            tabPanel("DeltaSASA",
              DT::dataTableOutput("popsDeltaSASAResidue"),
              tags$hr(),
              downloadButton('downloadResidueDeltaSASA', 'Download Residue DeltaSASA')
            ),
            tabPanel("Isolated Chains",
              DT::dataTableOutput("popsIsoSASAResidue"),
              tags$hr(),
              downloadButton('downloadResidueIsoSASA', 'Download Isolated-Chains Residue SASA')
            )
          )
        ),
        tabPanel("Chain",
          tabsetPanel(
            tabPanel("Input Structure",
              DT::dataTableOutput("popsSASAChain"),
              tags$hr(),
              downloadButton('downloadChainSASA', 'Download Chain SASA')
            ),
            tabPanel("DeltaSASA",
              DT::dataTableOutput("popsDeltaSASAChain"),
              tags$hr(),
              downloadButton('downloadChainDeltaSASA', 'Download Chain DeltaSASA')
            ),
            tabPanel("Isolated Chains",
              DT::dataTableOutput("popsIsoSASAChain"),
              tags$hr(),
              downloadButton('downloadChainIsoSASA', 'Download Isolated-Chains Chain SASA')
            )
          )
        ),
        tabPanel("Molecule",
          DT::dataTableOutput("popsMolecule"),
          tags$hr(),
          downloadButton('downloadMoleculeSASA', 'Download Molecule SASA')
        ),
        tabPanel("Usage",
		      h3("Method"),
          p("The POPScomp server invokes the POPS program to compute the
		        Solvent Accessible Surface Area (SASA) of a given PDB structure.
			      This is the POPS functionality, which has been currently implemented.
			      Soon to come: The POPSCOMP functionality for protein or RNA/DNA complexes,
			      where the POPScomp server creates internally
			      all pair combinations of chains to compute the buried SASA upon complexation.
			      Details of those functionalities are explained in the published papers
			      on POPS and POPSCOMP, see 'About' tab for the list of publications."
          ),
		      h3("Run"),
          p("SASA tables are initialised without any values; therefore, before 'run POPScomp' execution,
            the user sees only the table header and below the notice 'Showing 0 to 0 of 0 entries'.
            After selecting a PDB file and pressing 'run POPScomp', the server runs
            the POPS program on the selected PDB file. The output is SASA tables,
            which are automatically loaded into the respective tabs.
            Because that computation is a shell command, the success of the computation is
            returned as exit code and shown below the 'run POPScomp' button.
            It should normally read 'Exit code: 0', otherwise consult the 'Exit Codes' tab."
          ),
		      h3("Results"),
		      p("The SASA result tabs are 'Atom', 'Residue', 'Chain' and 'Molecule'.
            Except for 'Molecule', they all contain a second layer of tabs to accommodate
            the POPSCOMP functionality, as follows.
            'Input Structure': SASA values of the PDB structure as input.
            'DeltaSASA': The SASA difference between isolated chains and chain pair complexes.
            'Isolated Chains': SASA values of isolated chains.
            Only structures containing multiple chains will yield values for
            'DeltaSASA' and 'Isolated Chains' tabs."
		      ),
		      p("Results will be kept for one day on the server. Please use the 'Download ...' buttons
		        under the tables to save their content in 'csv' format. The 'Download All Results'
		        button on the side panel returns the zipped content of the output directory."
		      ),
		      h3("Help"),
          p("In case the server does not work as expected or server-related issues
		        need clarification, please email the maintainers:
			      Jens Kleinjung (jens@jkleinj.eu) and
            Franca Fraternali (franca.fraternali@kcl.ac.uk).
            For software and output errors, feature suggestions and similar topics,
            please add an entry to the ",
            a("Issues tab on the POPScomp GitHub page", href="https://github.com/Fraternalilab/POPScomp/issues"), "."
          )
        ),
        tabPanel("About",
			    h3("Servers"),
			    h4("Software"),
          p("This is version 3.0.6 of the", a("POPScomp server", href="http://popscom.org:3838"), "."),
			    p("Source code and detailed information will be released in 2020
			        on Fraternali Lab's GitHub page as ",
			      a("POPScomp repository", href="https://github.com/Fraternalilab/POPScomp"), "."
			    ),
			    p("The POPScomp server is based on two software packages:"),
          p("1. POPSC: A GNU Autotools package of the POPS C program."),
          p("2. POPSR: An R package containing this Shiny server
             to interface the POPS program and to provide the POPSCOMP functionality."
          ),
          p("Since April 2019, POPSC and POPSR are being co-developed."),
			    p("The legacy codes of POPS and POPSCOMP are available as repositories
			      'POPSlegacy' and 'POPSCOMPlegacy' on ",
			        a("Fraternali Lab's GitHub page", href="https://github.com/Fraternalilab"), "."
			    ),
			    h4("EBI PDBe-KB"),
          p("POPScomp is part of the ",
              a("FunPDBe resources", href="https://www.ebi.ac.uk/pdbe/funpdbe/deposition"), "."
          ),
			    h3("References"),
			    p("Users publishing results obtained with the program and
			        its applications should acknowledge its use by citation."),
			    h4("Implicit solvent"),
			    p("Fraternali, F. and van Gunsteren, W.F.
			        An efficient mean solvation force model for use in
			        molecular dynamics simulations of proteins in aqueous solution.
			        Journal of Molecular Biology 256 (1996) 939-948.",
			      a("DOI", href="https://dx.doi.org/10.1016%2Fj.sbi.2014.04.003"),
			      a("Pubmed", href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4045398/")
			    ),
			    h4("POPS method"),
			    p("Fraternali, F. and Cavallo, L.
			        Parameter optimized surfaces (POPS): analysis of key interactions
			        and conformational changes in the ribosome.
			        Nucleic Acids Research 30 (2002) 2950-2960.",
			      a("DOI", href="https://dx.doi.org/10.1093%2Fnar%2Fgkf373"),
			      a("Pubmed", href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC117037/")
			    ),
			    h4("POPS server"),
			    p("Cavallo, L., Kleinjung, J. and Fraternali, F.
			        POPS: A fast algorithm for solvent accessible surface areas
			        at atomic and residue level.
			        Nucleic Acids Research 31 (2003) 3364-3366.",
			      a("DOI", href="https://dx.doi.org/10.1093%2Fnar%2Fgkg601"),
			      a("Pubmed", href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC169007/")
			    ),
			    h4("POPSCOMP server"),
			    p("Kleinjung, J. and Fraternali, F.
			        POPSCOMP: an automated interaction analysis of biomolecular complexes.
			        Nucleic Acids Research 33 (2005) W342-W346.",
			      a("DOI", href="https://dx.doi.org/10.1093%2Fnar%2Fgki369"),
			      a("Pubmed", href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1160130/")
			    ),
			    h3("License"),
			    p("Usage of the software and server is free under the
			        GNU General Public License v3.0."
			    ),
			    h4("Copyright Holders, Authors and Maintainers"),
			    p("2002-2020 Franca Fraternali (author, maintainer)"),
			    p("2008-2020 Jens Kleinjung (author, maintainer)"),
			    h4("Contributors"),
			    p("2002 Kuang Lin and Valerie Hindie (translation to C)"),
			    p("2002 Luigi Cavallo (parametrisation)")
        ),
        tabPanel("Exit Codes",
          h3("Shiny exit codes"),
          p("* No PDB source input! - Enter PDB identifier or upload PDB file from local file system
            at the top of the side panel."
          ),
          p("* Two PDB sources input! - Only one PDB source is accepted per computation. Refresh the
            browser page and either speficy a PDB identifier or upload a PDB file, not both."
          ),
          h3("Shell command exit codes"),
          p("* 0 - Success"),
          p("* 1 - Catchall for general errors"),
          p("* 2 - Misuse of shell builtins (according to Bash documentation"),
          p("* 126 - Command invoked cannot execute"),
          p("* 127 - Command not found"),
          p("* 128 - Invalid argument to exit"),
          p("* 128+n - Fatal error signal 'n'"),
          p("* 130 - Script terminated by Control-C"),
          p("* 255* - Exit status out of range"),
          h3("Troubleshooting Errors"),
          h4("Exit code: 1 AND Error: Cannot open the connection"),
          p("The PDB file could not be read, most possibly because something went wrong during up/down-loading.
            If you used the 'Enter PDB entry' field, check your internet connection."
          )
        )
      )
    )
  )
)

#_______________________________________________________________________________
# server routines
server <- function(input, output) {

  ## random output paths
  rndString = as.character(digest(format(Sys.time(), "%H:%M:%OS3")))
  mainDir = "/tmp"
  subDir = paste0("POPScomp_", rndString)
  outDir = paste(mainDir, subDir, sep = "/")

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
  ##     to which the specified PDB file will be up/down-loaded.
  ## - For uploaded files: The 'fileInput' function returns the object 'input$file1',
  ##     a list of four elements, of which the fourth element contains the
  ##     path to the temporary file.
  ## - POPS will be run on the PDB file and the output will be zipped.
  output$nil <- eventReactive(input$popscomp, {
    ## to proceed, we require one PDB identifier or uploaded PDB file
    validate(need(((input$pdbentry != "") || (! is.null(input$PDBfile))),
          message = "No PDB source input!"))
    ## to proceed, we need one unspecified PDB identifier
    validate(need(((input$pdbentry == "") || (is.null(input$PDBfile))),
          message = "Two PDB sources input!"))

    #mainDir = "/tmp"
    #rndString = as.character(digest(format(Sys.time(), "%H:%M:%OS3")))
    #subDir = paste0("POPScomp_", rndString)
    #outDir = paste(mainDir, subDir, sep = "/")
    dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
    setwd(file.path(mainDir, subDir))

    ## download (PDB database) or upload (local file system) the PDB structure
    if (input$pdbentry != "") {
      ## get.pdb downloads the PDB structure from the database
      get.pdb(input$pdbentry, path = outDir);
      inputPDB = paste(input$pdbentry, "pdb", sep = ".")
    } else {
      ## move uploaded PDB file from its temporary directory to output directory
      system(paste("mv ", input$PDBfile[[4]], " ", outDir, "/",
                   input$PDBfile[[1]],  sep = ""))
      inputPDB = input$PDBfile[[1]]
    }

    ## run POPS as system command
    if (input$popsmode == "atomistic") {
      command = paste("/usr/local/bin/pops --outDirName", outDir,
                      "--rout --atomOut --residueOut --chainOut",
                      "--rProbe", input$rprobe, "--pdb", inputPDB, "1> POPScomp.o 2> POPScomp.e");
    } else if (input$popsmode == "coarse") {
      command = paste("/usr/local/bin/pops --outDirName", outDir,
                      "--rout --coarse --chainOut --residueOut --chainOut",
                      "--rProbe", input$rprobe, "--pdb", inputPDB, "1> POPScomp.o 2> POPScomp.e");
    }
    system_status = system(command)
    ## zip output directory for potential All-Result download
    zip(paste0(mainDir, "/", subDir, ".zip"), outDir)
    ## return exit code of POPS command
    paste("Exit code:", system_status)
  })

  ## o5.1.1 atom SASA
  ## empty dataframe with column names
  ## that will show up as empty table before POPS has finished
  atom_sasa_null.df = data.frame(
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
  write.table(atom_sasa_null.df, file = "id.rpopsAtom")
  ## reactive data: update output when file content changes
  atomSASAOutput = reactiveValues(highlight = NULL, data = NULL)
  atomSASAOutputData = reactiveFileReader(2000, NULL, "id.rpopsAtom",
                                      read.table, header = TRUE)
  ## render output data as table
  output$popsSASAAtom = DT::renderDataTable({
    atomSASAOutput$data = atomSASAOutputData()
  })

  ## o5.1.2 atom DeltaSASA
  atom_deltasasa_null.df = data.frame(
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
  write.table(atom_deltasasa_null.df, file = "dSASA.rpopsAtom")
  atomDeltaSASAOutput = reactiveValues(highlight = NULL, data = NULL)
  atomDeltaSASAOutputData = reactiveFileReader(2000, NULL, "dSASA.rpopsAtom",
                                      read.table, header = TRUE)
  output$popsDeltaSASAAtom = DT::renderDataTable({
    atomDeltaSASAOutput$data = atomDeltaSASAOutputData()
  })

  ## o5.1.3 atom isolated-chain SASA
  atom_isosasa_null.df = data.frame(
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
  write.table(atom_isosasa_null.df, file = "isoSASA.rpopsAtom")
  atomIsoSASAOutput = reactiveValues(highlight = NULL, data = NULL)
  atomIsoSASAOutputData = reactiveFileReader(2000, NULL, "isoSASA.rpopsAtom",
                                         read.table, header = TRUE)
  output$popsIsoSASAAtom = DT::renderDataTable({
    atomIsoSASAOutput$data = atomIsoSASAOutputData()
  })

  ## o5.2.1 residue SASA
  residue_sasa_null.df = data.frame(
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
  write.table(residue_sasa_null.df, file = "id.rpopsResidue")
  residueSASAOutput = reactiveValues(highlight = NULL, data = NULL)
  residueSASAOutputData = reactiveFileReader(2000, NULL, "id.rpopsResidue",
                                         read.table, header = TRUE)
  output$popsSASAResidue = DT::renderDataTable({
    residueSASAOutput$data = residueSASAOutputData()
  })

  ## o5.2.2 residue DeltaSASA
  residue_deltasasa_null.df = data.frame(
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
  write.table(residue_deltasasa_null.df, file = "dSASA.rpopsResidue")
  residueDeltaSASAOutput = reactiveValues(highlight = NULL, data = NULL)
  residueDeltaSASAOutputData = reactiveFileReader(2000, NULL, "dSASA.rpopsResidue",
                                         read.table, header = TRUE)
  output$popsDeltaSASAResidue = DT::renderDataTable({
    residueDeltaSASAOutput$data = residueDeltaSASAOutputData()
  })

  ## o5.2.3 residue isolated-chain SASA
  residue_isosasa_null.df = data.frame(
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
  write.table(residue_isosasa_null.df, file = "isoSASA.rpopsResidue")
  residueIsoSASAOutput = reactiveValues(highlight = NULL, data = NULL)
  residueIsoSASAOutputData = reactiveFileReader(2000, NULL, "isoSASA.rpopsResidue",
                                      read.table, header = TRUE)
  output$popsIsoSASAResidue = DT::renderDataTable({
    residueIsoSASAOutput$data = residueIsoSASAOutputData()
  })

  ## o5.3.1 chain SASA
  chain_sasa_null.df = data.frame(
    ChainNr = integer(),
    ChainNe = character()
  )
  write.table(chain_sasa_null.df, file = "id.rpopsChain")
  chainSASAOutput = reactiveValues(highlight = NULL, data = NULL)
  chainSASAOutputData = reactiveFileReader(2000, NULL, "id.rpopsChain",
                                       read.table, header = TRUE)
  output$popsSASAChain = DT::renderDataTable({
    chainSASAOutput$data = chainSASAOutputData()
  })

  ## o5.3.2 chain DeltaSASA
  chain_deltasasa_null.df = data.frame(
    ChainNr = integer(),
    ChainNe = character()
  )
  write.table(chain_deltasasa_null.df, file = "dSASA.rpopsChain")
  chainDeltaSASAOutput = reactiveValues(highlight = NULL, data = NULL)
  chainDeltaSASAOutputData = reactiveFileReader(2000, NULL, "dSASA.rpopsChain",
                                       read.table, header = TRUE)
  output$popsDeltaSASAChain = DT::renderDataTable({
    chainDeltaSASAOutput$data = chainDeltaSASAOutputData()
  })

  ## o5.3.3 chain isolated-chain SASA
  chain_isosasa_null.df = data.frame(
                    ChainNr = integer(),
                    ChainNe = character()
  )
  write.table(chain_isosasa_null.df, file = "isoSASA.rpopsChain")
  chainIsoSASAOutput = reactiveValues(highlight = NULL, data = NULL)
  chainIsoSASAOutputData = reactiveFileReader(2000, NULL, "isoSASA.rpopsChain",
                                      read.table, header = TRUE)
  output$popsIsoSASAChain = DT::renderDataTable({
    chainIsoSASAOutput$data = chainIsoSASAOutputData()
  })

  ## o5.4 molecule
  molecule_sasa_null.df = data.frame(
                      Phob = double(),
                      Phil = double(),
                      Total = double()
  )
  write.table(molecule_sasa_null.df, file = "id.rpopsMolecule")
  moleculeSAASOutput = reactiveValues(highlight = NULL, data = NULL)
  moleculeSASAOutputData = reactiveFileReader(2000, NULL, "id.rpopsMolecule",
                                      read.table, header = TRUE)
  output$popsSASAMolecule = DT::renderDataTable({
    moleculeSASAOutput$data = moleculeSASAOutputData()
  })

  ## d1 function removed
  ## d2 download all results
  output$downloadAllResults <- downloadHandler(
    filename = function() {
      paste0(mainDir, "/", subDir, ".zip")
    },
    content = function(file) {
      file.copy(paste0(mainDir, "/", subDir, ".zip"), file)
    },
    contentType = "application/zip"
  )

  ## d3.1 download atom SASA
  output$downloadAtomSASA <- downloadHandler(
    filename = function() {
      paste('atomSASA_', format(Sys.time(), "%Y%m%d%H%M"), '.csv', sep = '')
    },
    content = function(fname) {
      write.csv(atomSASAOutputData(), fname)
    }
  )

  ## d3.2 download atom DeltaSASA
  output$downloadAtomDeltaSASA <- downloadHandler(
    filename = function() {
      paste('atomDeltaSASA_', format(Sys.time(), "%Y%m%d%H%M"), '.csv', sep = '')
    },
    content = function(fname) {
      write.csv(atomDeltaSASAOutputData(), fname)
    }
  )

  ## d3.3 download atom isolated-chains SASA
  output$downloadAtomIsoSASA <- downloadHandler(
    filename = function() {
      paste('atomIsoSASA_', format(Sys.time(), "%Y%m%d%H%M"), '.csv', sep = '')
    },
    content = function(fname) {
      write.csv(atomIsoSASAOutputData(), fname)
    }
  )

  ## d4.1 download residue SASA
  output$downloadResidueSASA <- downloadHandler(
    filename = function() {
      paste('residueSASA_', format(Sys.time(), "%Y%m%d%H%M"), '.csv', sep = '')
    },
    content = function(fname) {
      write.csv(residueSASAOutputData(), fname)
    }
  )

  ## d4.2 download residue DeltaSASA
  output$downloadResidueDeltaSASA <- downloadHandler(
    filename = function() {
      paste('residueDeltaSASA_', format(Sys.time(), "%Y%m%d%H%M"), '.csv', sep = '')
    },
    content = function(fname) {
      write.csv(residueDeltaSASAOutputData(), fname)
    }
  )

  ## d4.1 download residue isolated-chains SASA
  output$downloadResidueIsoSASA <- downloadHandler(
    filename = function() {
      paste('residueIsoSASA_', format(Sys.time(), "%Y%m%d%H%M"), '.csv', sep = '')
    },
    content = function(fname) {
      write.csv(residueIsoSASAOutputData(), fname)
    }
  )

  ## d5.1 download chain SASA
  output$downloadChainSASA <- downloadHandler(
    filename = function() {
      paste('chainSASA_', format(Sys.time(), "%Y%m%d%H%M"), '.csv', sep = '')
    },
    content = function(fname) {
      write.csv(chainSASAOutputData(), fname)
    }
  )

  ## d5.2 download chain DeltaSASA
  output$downloadChainDeltaSASA <- downloadHandler(
    filename = function() {
      paste('chainDeltaSASA_', format(Sys.time(), "%Y%m%d%H%M"), '.csv', sep = '')
    },
    content = function(fname) {
      write.csv(chainDeltaSASAOutputData(), fname)
    }
  )

  ## d5.1 download chain isolated-chains SASA
  output$downloadChainIsoSASA <- downloadHandler(
    filename = function() {
      paste('chainIsoSASA_', format(Sys.time(), "%Y%m%d%H%M"), '.csv', sep = '')
    },
    content = function(fname) {
      write.csv(chainIsoSASAOutputData(), fname)
    }
  )

  ## d6 download molecule SASA
  output$downloadMoleculeSASA <- downloadHandler(
    filename = function() {
      paste('moleculeSASA_', format(Sys.time(), "%Y%m%d%H%M"), '.csv', sep = '')
    },
    content = function(fname) {
      write.csv(moleculeSASAOutputData(), fname)
    }
  )
}

#_______________________________________________________________________________
# run the Shiny app
shinyApp(ui, server)

#===============================================================================
