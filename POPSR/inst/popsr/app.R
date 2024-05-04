#===============================================================================
# Shiny application as interface of POPScomp
# Here the C program POPS is called directly only on the original input file
#   (tagged 'iso'). All additional POPS(comp) calls/computations
#   are related to processing protein complexes,
#   and those are coded in .../POPScomp/R/popscomp.R.
# For single-chain proteins, "popscomp.R" returns without additional computations.
#
# (C) 2019 -2023 Jens Kleinjung and Franca Fraternali
#===============================================================================

library(shiny)
library(bio3d)
library(DT)
library(digest)
library(shinysky)
library(POPSR)

#_______________________________________________________________________________
## load 'Readme' text
readme = readRDS("readme.rds")

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
      fileInput(inputId = "PDBfile",
                label = "OR upload PDB file",
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

      ## info1 run ID
      textOutput("runid"),
      tags$hr(),

      busyIndicator(wait = 1000),

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
          tabsetPanel(
            tabPanel("Input Structure",
              DT::dataTableOutput("popsSASAMolecule"),
              tags$hr(),
              downloadButton('downloadMoleculeSASA', 'Download Molecule SASA')
            ),
            tabPanel("DeltaSASA",
              DT::dataTableOutput("popsDeltaSASAMolecule"),
              tags$hr(),
              downloadButton('downloadMoleculeDeltaSASA', 'Download Molecule DeltaSASA')
            )
          )
        ),
        tabPanel("Usage",
		      h3("Method"),
          p("The POPScomp server invokes the POPS program to compute the
		        Solvent Accessible Surface Area (SASA) of a given PDB structure.
			      For protein or RNA/DNA complexes, the POPScomp server creates internally
			      all pair combinations of chains to compute the buried SASA upon complexation.
			      Details of those functionalities are explained in the published papers
			      on implicit solvent, POPS and POPSCOMP; see the 'About' tab for the list of publications."
          ),
		      h3("Run"),
          p("SASA tables are initialised without any values; therefore, before 'run POPScomp' execution,
            the user sees only the table header and below the notice 'Showing 0 to 0 of 0 entries'.
            After selecting a PDB identifier or uploading a PDB file and pressing 'run POPScomp', the server runs
            the POPS program on the PDB input. The output are SASA tables,
            which are automatically loaded into the respective tabs on the POPScomp interface.
            The success of the computation is returned as exit code and shown below
	        the 'run POPScomp' button: 'Exit code: 0' means success and that is what you
	        should expect to see, otherwise consult the 'Exit Codes' tab.
            The 'run ID' identifier is a random string that is updated upon changes in the
            input parameters. The identifier is used in the ouput path for the results,
            Therefore, when re-running POPScomp on the same PDB input, changed parameters will always yield
            a separate output, whereas unchanged settings (identical random string) will overwrite the previous results."
          ),
		      h3("Results"),
		      p("The SASA result tabs are 'Atom', 'Residue', 'Chain' and 'Molecule'.
                Those tabs contain a second layer of tabs to accommodate POPSCOMP's complex analysis, as follows."),
              p("'Input Structure': SASA values of the input PDB structure."),
              p("'DeltaSASA': The SASA difference between isolated chains and chain pair complexes.
					The DeltaSASA values correspond to the buried surface area upon complexation between the two chains
                    in the given chain pair."
			  ),
              p("'Isolated Chains': SASA values of isolated chains. These values form the basis of the DeltaSASA data."),
              p("Only structures containing multiple chains will yield values for 'DeltaSASA' and 'Isolated Chains' tabs."),
		      p("Results will be stored on the server until midnight GMT time and then automatically removed.
		        Please use the 'Download ...' buttons under the tables to save your results in 'csv' format.
		        The 'Download All Results' button on the side panel returns the zipped content of the entire output directory,
		        i.e. all results produced for a given POPScomp job."
		      ),
			h3("Output Columns"),
			h4("ATOM SASAs"),
			p("AtomNr : atom number in molecular coordinate file"),
			p("AtomNe : atom name in molecular coordinate file"),
			p("ResidNe : residue name in molecular coordinate file"),
			p("Chain : chain name in molecular coordinate file ('-' if unspecified)"),
			p("ResidNr : residue number in molecular coordinate file"),
			p("SASA.A.2 : solvent accessible surace area in Angstrom^2 units"),
			p("Q.SASA. : quotient of SASA and Surf (below), i.e. the fraction of SASA"),
			p("N.overl. : number of overlaps with atom neighbours"),
			p("AtomTp : atom type code (GROMOS van der Waals atom type)"),
			p("AtomGp : atom group code: positive=1, negative=2, polar=3, aromatic=4, aliphatic=5"),
			p("Surf/A^2 : surface area of isolated atom"),

			h4("RESIDUE SASAs"),
			p("ResidNe : residue name in molecular coordinate file"),
			p("Chain : chain name in molecular coordinate file ('-' if unspecified)"),
			p("ResidNr : residue number in molecular coordinate file"),
			p("Phob.A.2 : hydrophobic solvent accessible surace area in Angstrom^2 units"),
			p("Phil.A.2 : hydrophilic solvent accessible surace area in Angstrom^2 units"),
			p("SASA.A.2 : total solvent accessible surace area in Angstrom^2 units"),
			p("Q.SASA. : quotient of SASA and Surf (below), i.e. the fraction of SASA"),
			p("N.overl. : number of overlaps with residue neighbours"),
			p("Surf.A.2 : surface area of isolated residue"),

			h4("CHAIN SASAs"),
			p("Chain : chain number"),
			p("Id : chain name in molecular coordinate file ('-' if unspecified)"),
			p("AtomRange : range of atom numbers in chain"),
			p("ResidRange : range of residue numbers in chain"),
			p("Phob.A.2 : hydrophobic solvent accessible surace area in Angstrom^2 units"),
			p("Phil.A.2 : hydrophilic solvent accessible surace area in Angstrom^2 units"),
			p("SASA.A.2 : total solvent accessible surace area in Angstrom^2 units"),

			h4("MOLECULE SASAs"),
			p("Phob.A.2 : hydrophobic solvent accessible surace area in Angstrom^2 units"),
			p("Phil.A.2 : hydrophilic solvent accessible surace area in Angstrom^2 units"),
			p("SASA.A.2 : total solvent accessible surace area in Angstrom^2 units"),

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
			h3("Shiny App"),
			p("This is version 3.3 of the POPScomp Shiny App."),
			p("For detailed information about the software visit Fraternali Lab's ",
			  a("POPScomp GitHub repository", href="https://github.com/Fraternalilab/POPScomp"),
			  "; the Wiki pages contain detailed installation and usage instructions."
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
			p("Kleinjung, J. and Fraternali, F.
			Design and Application of Implicit Solvent Models
			in Biomolecular Simulations.
                            Current Opinion in Structural Biology 25 (2014) 126-134.",
			  a("DOI", href="http://dx.doi.org/10.1016/j.sbi.2014.04.003"),
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
			h3("License and Copyright"),
			p("Usage of the software and server is free under the
			    GNU General Public License v3.0."
			),
			h4("Copyright Holders, Authors and Maintainers"),
			p("2002-2023 Franca Fraternali (author, maintainer)"),
			p("2008-2023 Jens Kleinjung (author, maintainer)"),
			h4("Contributors"),
			p("2002 Kuang Lin and Valerie Hindie (translation to C)"),
			p("2002 Luigi Cavallo (parametrisation)")
        ),
        tabPanel("Exit Codes",
			h3("Overview"),
			p("POPScomp uses a combination of *Shell* (system) calls and R *Shiny* routines.
				Therefore, the return value shown as exit code may come from *Shell* or *Shiny*.
				A successful run will return 'Exit code: 0'. Any error will return an exit code
				different from '0'. A commented list of exit codes is given below together with
				troubleshooting tips. In case you get stuck, please contact the maintainers."),
          h3("Shell command exit codes"),
          p("* 0 - Success"),
          p("* 1 - Catchall for general errors"),
          p("* 2 - Misuse of shell builtins (according to Bash documentation)"),
          p("* 126 - Command invoked cannot execute"),
          p("* 127 - Command not found"),
          p("* 128 - Invalid argument to exit"),
          p("* 128+n - Fatal error signal 'n'"),
          p("* 130 - Script terminated by Control-C"),
          p("* 255* - Exit status out of range"),
          h3("Shiny exit codes"),
          p("* No PDB source input! - Enter PDB identifier or upload PDB file from local file system
            at the top of the side panel."
          ),
          p("* Two PDB sources input! - Only one PDB source is accepted per computation. Refresh the
            browser page and either specify a PDB identifier or upload a PDB file, not both."
          ),
          h3("Troubleshooting Errors"),
          h4("Exit code: 1 AND Error: Cannot open the connection"),
          p("The PDB file could not be read, most possibly because something went wrong during up/down-loading.
            If you used the 'Enter PDB entry' field, check your internet connection."
          )
        ),
			tabPanel("Readme",
			    h3("Readme (C code)"),
			    textOutput("readme")
			  )
      )
    )
  )
)

#_______________________________________________________________________________
# server routines
server <- function(input, output) {

  ## pops binary in Shiny installation via Dockerfile
  pops_shiny = c("/build/install/usr/local/bin/pops")

  ## initialisation output directory
  ## Shiny needs empty output files to be present here,
  ##   otherwise it shows error messages on the start page before the POPScomp run.
  ## This will be overwritten/replaced by the input$popscomp function,
  ##   which will copy over the empty initialisation files and overwrite them where applicable.
  mainDir = "/tmp"
  subDir = "POPScomp_init"
  outDir = paste(mainDir, subDir, sep = "/")
  dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
  setwd(file.path(mainDir, subDir))

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

  ## info1 run identifier
  ## creation of the run identifier is reactive to any input parameter change
  ## the run identifier is a digest of the system time
  runid_string <- eventReactive({input$pdbentry
                                 input$PDBfile
                                 input$popsmode
                                 input$rprobe
                                 },{
    as.character(digest(format(Sys.time(), "%H%M%OS3")))
  })
  output$runid <- renderText({
    paste("run ID: ", runid_string(), sep = '')
  })

  ## o4 download PDB entry or upload input file
  ## run POPS on specified PDB file
  ## Comments:
  ## - 'pops' binary defined as 'pops_shiny' at beginning of server section.
  ## - The App will set its own working directory to a temporary directory,
  ##     to which the specified PDB file will be up/down-loaded.
  ## - For uploaded files: The 'fileInput' function returns the object 'input$file1',
  ##     a list of four elements, of which the fourth element contains the
  ##     path to the temporary file.
  ## - POPS will be run on the PDB file and the output will be zipped.
  output$nil <- eventReactive(input$popscomp, {
    ## create a new output directory for each POPScomp run
    subDir = paste0("POPScomp_", runid_string())
    outDir = paste(mainDir, subDir, sep = "/")
    dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
    setwd(file.path(mainDir, subDir))

    ## Copy initialisation files to the new output directory,
    ##   otherwise error messages will appear in non-complex structure results
    ##   for isolated chains and difference SASAs.
    initDir = "POPScomp_init"
    cpInit = paste0("cp ", mainDir, "/", initDir, "/* .")
    system(cpInit, wait = TRUE)

    ## to proceed, we require one PDB identifier or uploaded PDB file
    validate(need(((input$pdbentry != "") || (! is.null(input$PDBfile))),
          message = "No PDB source input!"))
    ## to proceed, we need one unspecified PDB identifier
    validate(need(((input$pdbentry == "") || (is.null(input$PDBfile))),
          message = "Two PDB sources input!"))

    ## download (PDB database) or upload (local file system) the PDB structure
    if (input$pdbentry != "") {
      ## get.pdb downloads the PDB structure from the database
      get.pdb(input$pdbentry, path = outDir);
      inputPDB = paste(input$pdbentry, "pdb", sep = ".")
    } else {
      ## move uploaded PDB file from its temporary directory to output directory
      system(paste("mv ", input$PDBfile[[4]], " ", outDir, "/",
                   input$PDBfile[[1]],  sep = ""), wait = TRUE)
      inputPDB = input$PDBfile[[1]]
    }

    ## run POPS as system command
    if (input$popsmode == "atomistic") {
      command = paste(pops_shiny, "--outDirName", outDir,
                      "--rout --atomOut --residueOut --chainOut",
                      "--rProbe", input$rprobe, "--pdb", inputPDB, "1> POPScomp.o 2> POPScomp.e");
    } else if (input$popsmode == "coarse") {
      command = paste(pops_shiny, "--outDirName", outDir,
                      "--rout --coarse --chainOut --residueOut --chainOut",
                      "--rProbe", input$rprobe, "--pdb", inputPDB, "1> POPScomp.o 2> POPScomp.e");
    }
    system_status = system(command, wait = TRUE)

    ## run POPScomp
    popscompR(inputPDB, outDir)
    ## zip output directory for potential All-Result download
    zip(paste0("POPScomp_", runid_string()), outDir)
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
                        SASA.A.2 = double(),
                        Q.SASA. = double(),
                        N.overl. = integer(),
                        AtomTp = integer(),
                        AtomGp = integer(),
                        Surf.A.2 = double()
  )
  write.table(atom_sasa_null.df, file = paste(outDir, "id.rpopsAtom", sep = '/'))
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
                            D_SASA.A.2 = double(),
                            AtomTp = integer(),
                            AtomGp = integer()
  )
  write.table(atom_deltasasa_null.df, file = paste(outDir, "deltaSASA.rpopsAtom", sep = '/'))
  atomDeltaSASAOutput = reactiveValues(highlight = NULL, data = NULL)
  atomDeltaSASAOutputData = reactiveFileReader(2000, NULL, "deltaSASA.rpopsAtom",
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
                          SASA.A.2 = double(),
                          Q.SASA. = double(),
                          N.overl. = integer(),
                          AtomTp = integer(),
                          AtomGp = integer(),
                          Surf.A.2 = double()
  )
  write.table(atom_isosasa_null.df, file = paste(outDir, "isoSASA.rpopsAtom", sep = '/'))
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
                          Phob.A.2 = double(),
                          Phil.A.2 = double(),
                          SASA.A.2 = double(),
                          Q.SASA. = double(),
                          N.overl. = integer(),
                          Surf.A.2 = double()
  )
  write.table(residue_sasa_null.df, file = paste(outDir, "id.rpopsResidue", sep = '/'))
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
                                D_Phob.A.2 = double(),
                                D_Phil.A.2 = double(),
                                D_SASA.A.2 = double()
  )
  write.table(residue_deltasasa_null.df, file = paste(outDir, "deltaSASA.rpopsResidue", sep = '/'))
  residueDeltaSASAOutput = reactiveValues(highlight = NULL, data = NULL)
  residueDeltaSASAOutputData = reactiveFileReader(2000, NULL, "deltaSASA.rpopsResidue",
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
                              Phob.A.2 = double(),
                              Phil.A.2 = double(),
                              SASA.A.2 = double(),
                              Q.SASA. = double(),
                              N.overl. = integer(),
                              Surf.A.2 = double()
  )
  write.table(residue_isosasa_null.df, file = paste(outDir, "isoSASA.rpopsResidue", sep = '/'))
  residueIsoSASAOutput = reactiveValues(highlight = NULL, data = NULL)
  residueIsoSASAOutputData = reactiveFileReader(2000, NULL, "isoSASA.rpopsResidue",
                                      read.table, header = TRUE)
  output$popsIsoSASAResidue = DT::renderDataTable({
    residueIsoSASAOutput$data = residueIsoSASAOutputData()
  })

  ## o5.3.1 chain SASA
  chain_sasa_null.df = data.frame(
                          Chain = integer(),
                          Id = character(),
                          AtomRange = character(),
                          ResidRange = character(),
                          Phob.A.2 = double(),
                          Phil.A.2 = double(),
                          SASA.A.2 = double()
  )
  write.table(chain_sasa_null.df, file = paste(outDir, "id.rpopsChain", sep = '/'))
  chainSASAOutput = reactiveValues(highlight = NULL, data = NULL)
  chainSASAOutputData = reactiveFileReader(2000, NULL, "id.rpopsChain",
                                       read.table, header = TRUE)
  output$popsSASAChain = DT::renderDataTable({
    chainSASAOutput$data = chainSASAOutputData()
  })

  ## o5.3.2 chain DeltaSASA
  chain_deltasasa_null.df = data.frame(
                              Chain = integer(),
                              Id = character(),
                              AtomRange = character(),
                              ResidRange = character(),
                              D_Phob.A.2 = double(),
                              D_Phil.A.2 = double(),
                              D_SASA.A.2 = double()
  )
  write.table(chain_deltasasa_null.df, file = paste(outDir, "deltaSASA.rpopsChain", sep = '/'))
  chainDeltaSASAOutput = reactiveValues(highlight = NULL, data = NULL)
  chainDeltaSASAOutputData = reactiveFileReader(2000, NULL, "deltaSASA.rpopsChain",
                                       read.table, header = TRUE)
  output$popsDeltaSASAChain = DT::renderDataTable({
    chainDeltaSASAOutput$data = chainDeltaSASAOutputData()
  })

  ## o5.3.3 chain isolated-chain SASA
  chain_isosasa_null.df = data.frame(
                            Chain = integer(),
                            Id = character(),
                            AtomRange = character(),
                            ResidRange = character(),
                            Phob.A.2 = double(),
                            Phil.A.2 = double(),
                            SASA.A.2 = double()
  )
  write.table(chain_isosasa_null.df, file = paste(outDir, "isoSASA.rpopsChain", sep = '/'))
  chainIsoSASAOutput = reactiveValues(highlight = NULL, data = NULL)
  chainIsoSASAOutputData = reactiveFileReader(2000, NULL, "isoSASA.rpopsChain",
                                      read.table, header = TRUE)
  output$popsIsoSASAChain = DT::renderDataTable({
    chainIsoSASAOutput$data = chainIsoSASAOutputData()
  })

  ## o5.4.1 molecule SASA
  molecule_sasa_null.df = data.frame(
                            Phob.A.2 = double(),
                            Phil.A.2 = double(),
                            SASA.A.2 = double()
  )
  write.table(molecule_sasa_null.df, file = paste(outDir, "id.rpopsMolecule", sep = '/'))
  moleculeSASAOutput = reactiveValues(highlight = NULL, data = NULL)
  moleculeSASAOutputData = reactiveFileReader(2000, NULL, "id.rpopsMolecule",
                                      read.table, header = TRUE)
  output$popsSASAMolecule = DT::renderDataTable({
    moleculeSASAOutput$data = moleculeSASAOutputData()
  })

  ## o5.4.2 molecule DeltaSASA
  molecule_deltasasa_null.df = data.frame(
                                D_Phob.A.2 = double(),
                                D_Phil.A.2 = double(),
                                D_SASA.A.2 = double()
  )
  write.table(molecule_deltasasa_null.df, file = paste(outDir, "deltaSASA.rpopsMolecule", sep = '/'))
  moleculeDeltaSASAOutput = reactiveValues(highlight = NULL, data = NULL)
  moleculeDeltaSASAOutputData = reactiveFileReader(2000, NULL, "deltaSASA.rpopsMolecule",
                                                read.table, header = TRUE)
  output$popsDeltaSASAMolecule = DT::renderDataTable({
    moleculeDeltaSASAOutput$data = moleculeDeltaSASAOutputData()
  })

  ## d1 function removed
  ## d2 download all results
  output$downloadAllResults <- downloadHandler(
    filename = function() {
      paste0("POPScomp_", runid_string(),".zip")
    },
    content = function(file) {
      file.copy(paste0("POPScomp_", runid_string(),".zip"), file)
    },
    contentType = "application/zip"
  )

  ## d3.1 download atom SASA
  output$downloadAtomSASA <- downloadHandler(
    filename = function() {
      paste('atomSASA_', runid_string(), '.csv', sep = '')
    },
    content = function(fname) {
      write.csv(atomSASAOutputData(), fname)
    }
  )

  ## d3.2 download atom DeltaSASA
  output$downloadAtomDeltaSASA <- downloadHandler(
    filename = function() {
      paste('atomDeltaSASA_', runid_string(), '.csv', sep = '')
    },
    content = function(fname) {
      write.csv(atomDeltaSASAOutputData(), fname)
    }
  )

  ## d3.3 download atom isolated-chains SASA
  output$downloadAtomIsoSASA <- downloadHandler(
    filename = function() {
      paste('atomIsoSASA_', runid_string(), '.csv', sep = '')
    },
    content = function(fname) {
      write.csv(atomIsoSASAOutputData(), fname)
    }
  )

  ## d4.1 download residue SASA
  output$downloadResidueSASA <- downloadHandler(
    filename = function() {
      paste('residueSASA_', runid_string(), '.csv', sep = '')
    },
    content = function(fname) {
      write.csv(residueSASAOutputData(), fname)
    }
  )

  ## d4.2 download residue DeltaSASA
  output$downloadResidueDeltaSASA <- downloadHandler(
    filename = function() {
      paste('residueDeltaSASA_', runid_string(), '.csv', sep = '')
    },
    content = function(fname) {
      write.csv(residueDeltaSASAOutputData(), fname)
    }
  )

  ## d4.1 download residue isolated-chains SASA
  output$downloadResidueIsoSASA <- downloadHandler(
    filename = function() {
      paste('residueIsoSASA_', runid_string(), '.csv', sep = '')
    },
    content = function(fname) {
      write.csv(residueIsoSASAOutputData(), fname)
    }
  )

  ## d5.1 download chain SASA
  output$downloadChainSASA <- downloadHandler(
    filename = function() {
      paste('chainSASA_', runid_string(), '.csv', sep = '')
    },
    content = function(fname) {
      write.csv(chainSASAOutputData(), fname)
    }
  )

  ## d5.2 download chain DeltaSASA
  output$downloadChainDeltaSASA <- downloadHandler(
    filename = function() {
      paste('chainDeltaSASA_', runid_string(), '.csv', sep = '')
    },
    content = function(fname) {
      write.csv(chainDeltaSASAOutputData(), fname)
    }
  )

  ## d5.1 download chain isolated-chains SASA
  output$downloadChainIsoSASA <- downloadHandler(
    filename = function() {
      paste('chainIsoSASA_', runid_string(), '.csv', sep = '')
    },
    content = function(fname) {
      write.csv(chainIsoSASAOutputData(), fname)
    }
  )

  ## d6 download molecule SASA
  output$downloadMoleculeSASA <- downloadHandler(
    filename = function() {
      paste('moleculeSASA_', runid_string(), '.csv', sep = '')
    },
    content = function(fname) {
      write.csv(moleculeSASAOutputData(), fname)
    }
  )

  ## d7 download molecule DeltaSASA
  output$downloadMoleculeDeltaSASA <- downloadHandler(
    filename = function() {
      paste('moleculeDeltaSASA_', runid_string(), '.csv', sep = '')
    },
    content = function(fname) {
      write.csv(moleculeDeltaSASAOutputData(), fname)
    }
  )

  ## Readme panel text output (loaded as variable 'readme' at top of App)
  output$readme = renderText(readme)
}

#_______________________________________________________________________________
# run the Shiny app
shinyApp(ui, server)

#===============================================================================
