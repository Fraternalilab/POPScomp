#===============================================================================
# Shiny application as interface of POPScomp
# Here the C program POPS is called directly only on the original input file
#   (tagged 'iso'). All additional POPS(comp) calls/computations
#   are related to processing protein complexes,
#   and those are coded in .../POPScomp/R/popscomp.R.
# For single-chain proteins, "popscomp.R" returns without additional computations.
#
# (C) 2019-2026 Jens Kleinjung and Franca Fraternali
#===============================================================================

library(shiny)
library(bio3d)
library(DT)
library(digest)
library(shinysky)
library(POPSR)
library(markdown)

#_______________________________________________________________________________
## load 'Readme' text (optional: only used by the commented-out Readme tab)
readme = if (file.exists("readme.rds")) readRDS("readme.rds") else NULL

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
		      p("Please use the 'Download ...' buttons under the tables to save your results in 'csv' format.
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
          p("In case the program does not work as expected or server-related issues
		        need clarification, please email the maintainers:
			    Jens Kleinjung (jens@jkleinj.eu) and
                Franca Fraternali (f.fraternali@ucl.ac.uk).
            For software and output errors, feature suggestions and similar topics,
            please add an entry to the ",
            a("Issues tab on the POPScomp GitHub page", href="https://github.com/Fraternalilab/POPScomp/issues"), "."
          )
        ),
        tabPanel("About",
			h3("Shiny App"),
			p("This is version 3.5 of the POPScomp Shiny App."),
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
			p("2002-2026 Franca Fraternali (author, maintainer)"),
			p("2008-2026 Jens Kleinjung (author, maintainer)"),
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
        )
      )
    )
  )
)

#_______________________________________________________________________________
# server routines
server <- function(input, output, session) {

  ## Path of the POPS program: the POPS_BIN environment variable, else the
  ## program on the PATH, else the location in the Shiny Docker image.
  pops_bin = Sys.getenv("POPS_BIN", unset = "")
  if (! nzchar(pops_bin)) pops_bin = Sys.which("pops")
  if (! nzchar(pops_bin)) pops_bin = "/build/install/usr/local/bin/pops"

  ## Every session computes in its own directory and every file is addressed by
  ##   its absolute path. The working directory of the R process is never changed:
  ##   it is shared by all sessions, so one session would otherwise serve its
  ##   results into the browser of another.
  sessionDir = file.path(tempdir(), paste0("POPScomp_", digest(Sys.time())))
  dir.create(sessionDir, showWarnings = FALSE, recursive = TRUE)
  ## directory of the current run; before the first run there is none and the
  ##   result tables show their empty placeholders
  runDir = reactiveVal(sessionDir)
  ## run identifier of the last completed run, for the download file names
  runid_done = reactiveVal(NULL)
  ## remove the session directory when the browser tab is closed
  session$onSessionEnded(function() {
    unlink(sessionDir, recursive = TRUE)
  })

  ## Read one POPS output table. A table that does not exist yet, or that POPS
  ##   is in the middle of writing, yields the empty placeholder instead of an
  ##   error message in the user interface.
  readSasaTable = function(path, empty) {
    if (! file.exists(path)) return(empty)
    out = tryCatch({
        ## identifying columns are read as character: 'T' and 'F' are valid
        ##   chain identifiers and would otherwise become TRUE and FALSE
        header = scan(path, what = "", nlines = 1, quiet = TRUE)
        colClasses = rep(NA_character_, length(header))
        colClasses[header %in% c("AtomNe", "ResidNe", "Chain", "Id", "iCode",
                                 "AtomRange", "ResidRange")] = "character"
        read.table(path, header = TRUE, stringsAsFactors = FALSE,
                   colClasses = colClasses)
      }, error = function(e) empty)
    if (is.null(out) || nrow(out) == 0) return(empty)
    return(out)
  }

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
  ## - 'pops' binary located as 'pops_bin' at the beginning of the server section.
  ## - The App will set its own working directory to a temporary directory,
  ##     to which the specified PDB file will be up/down-loaded.
  ## - For uploaded files: The 'fileInput' function returns the object 'input$file1',
  ##     a list of four elements, of which the fourth element contains the
  ##     path to the temporary file.
  ## - POPS will be run on the PDB file and the output will be zipped.
  ## The run is an eventReactive so that it is evaluated when the result text is
  ##   rendered; 'validate' inside it reports its message in that text field.
  popscomp_run <- eventReactive(input$popscomp, {
    ## to proceed, we require one PDB identifier or uploaded PDB file
    ## (checked before anything is created, so that a failed attempt leaves the
    ##  results of the previous run in place)
    validate(need(((input$pdbentry != "") || (! is.null(input$PDBfile))),
          message = "No PDB source input!"))
    ## to proceed, we refuse 2 specified PDB inputs
    validate(need(((input$pdbentry == "") || (is.null(input$PDBfile))),
          message = "Two PDB sources input!"))
    ## a PDB identifier has four alphanumeric characters
    if (input$pdbentry != "") {
      validate(need(grepl("^[0-9A-Za-z]{4}$", input$pdbentry),
          message = "Not a PDB identifier: four letters or digits expected"))
    }

    ## create a new output directory for each POPScomp run
    runid = runid_string()
    outDir = file.path(sessionDir, paste0("POPScomp_", runid))
    dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

    ## download (PDB database) or upload (local file system) the PDB structure
    if (input$pdbentry != "") {
      ## get.pdb downloads the PDB structure from the database
      get.pdb(input$pdbentry, format = "cif", path = outDir)
      mmcifDownloadName = file.path(outDir, paste0(input$pdbentry, ".cif"))
      validate(need(file.exists(mmcifDownloadName),
          message = paste0("No structure ", input$pdbentry, " in the PDB")))
      inputPDB = paste0(input$pdbentry, ".pdb")
      command0 = paste("gemmi convert", shQuote(mmcifDownloadName),
                       shQuote(file.path(outDir, inputPDB)))
      system_status0 = system(command0)
      validate(need(system_status0 == 0,
          message = paste0("Conversion of ", input$pdbentry, " failed")))
    } else {
      ## copy the uploaded PDB file from its temporary directory to the output
      ##   directory; the name comes from the browser, so it is reduced to
      ##   harmless characters and the file is copied, not moved, so that the
      ##   same upload can be run more than once
      inputPDB = gsub("[^A-Za-z0-9._-]", "_", basename(input$PDBfile$name))
      validate(need(file.copy(input$PDBfile$datapath, file.path(outDir, inputPDB),
                              overwrite = TRUE),
          message = "Could not read the uploaded file"))
    }

    ## run POPS as system command
    coarse = (input$popsmode == "coarse")
    command = paste(shQuote(pops_bin), "--outDirName", shQuote(outDir),
                    "--rout --residueOut --chainOut --neighbourOut",
                    if (coarse) "--coarse" else "--atomOut",
                    "--rProbe", shQuote(as.character(input$rprobe)),
                    "--pdb", shQuote(file.path(outDir, inputPDB)),
                    "1>", shQuote(file.path(outDir, "POPScomp.o")),
                    "2>", shQuote(file.path(outDir, "POPScomp.e")))
    system_status = system(command, wait = TRUE)
    validate(need(system_status == 0,
        message = paste0("POPS declined this structure (exit code ", system_status,
                         "); see the Exit Codes tab")))

    ## run POPScomp on the chains and chain pairs of a complex
    popscompR(inputPDB, outDir, coarse = coarse)

    ## zip output directory for potential All-Result download;
    ##   the archive is written outside the directory it archives
    ## 'zip' is called with absolute paths and '-j' (no directory names), so that
    ##   the working directory of the process does not have to be changed
    zipFile = file.path(sessionDir, paste0("POPScomp_", runid, ".zip"))
    zip(zipFile, list.files(outDir, full.names = TRUE), flags = "-r9Xjq")

    ## publish the results of this run to the tables and the download buttons
    runDir(outDir)
    runid_done(runid)

    ## return exit code of POPS command
    paste("Exit code:", system_status)
  })

  output$nil <- renderText({
    popscomp_run()
  })

  ## o5.1.1 atom SASA
  ## empty dataframe with column names
  ## that will show up as empty table before POPS has finished
  atom_sasa_null.df = data.frame(
                        AtomNr = integer(),
                        AtomNe = character(),
                        ResidNe = character(),
                        Chain = character(),
                        ResidNr = integer(),
                        iCode = character(),
                        SASA.A.2 = double(),
                        Q.SASA. = double(),
                        N.overl. = integer(),
                        AtomTp = integer(),
                        AtomGp = integer(),
                        Surf.A.2 = double()
  )
  ## reactive data: update output when file content changes
  atomSASAOutput = reactiveValues(highlight = NULL, data = NULL)
  atomSASAOutputData = reactiveFileReader(2000, session,
                         function() file.path(runDir(), "id.rpopsAtom"),
                         readSasaTable, empty = atom_sasa_null.df)
  ## render output data as table
  output$popsSASAAtom = DT::renderDataTable({
    atomSASAOutput$data = atomSASAOutputData()
  })

  ## o5.1.2 atom DeltaSASA
  atom_deltasasa_null.df = data.frame(
                            AtomNr = integer(),
                            AtomNe = character(),
                            ResidNe = character(),
                            Chain = character(),
                            ResidNr = integer(),
                            iCode = character(),
                            D_SASA.A.2 = double(),
                            AtomTp = integer(),
                            AtomGp = integer()
  )
  atomDeltaSASAOutput = reactiveValues(highlight = NULL, data = NULL)
  atomDeltaSASAOutputData = reactiveFileReader(2000, session,
                         function() file.path(runDir(), "deltaSASA.rpopsAtom"),
                         readSasaTable, empty = atom_deltasasa_null.df)
  output$popsDeltaSASAAtom = DT::renderDataTable({
    atomDeltaSASAOutput$data = atomDeltaSASAOutputData()
  })

  ## o5.1.3 atom isolated-chain SASA
  atom_isosasa_null.df = data.frame(
                          AtomNr = integer(),
                          AtomNe = character(),
                          ResidNe = character(),
                          Chain = character(),
                          ResidNr = integer(),
                          iCode = character(),
                          SASA.A.2 = double(),
                          Q.SASA. = double(),
                          N.overl. = integer(),
                          AtomTp = integer(),
                          AtomGp = integer(),
                          Surf.A.2 = double()
  )
  atomIsoSASAOutput = reactiveValues(highlight = NULL, data = NULL)
  atomIsoSASAOutputData = reactiveFileReader(2000, session,
                         function() file.path(runDir(), "isoSASA.rpopsAtom"),
                         readSasaTable, empty = atom_isosasa_null.df)
  output$popsIsoSASAAtom = DT::renderDataTable({
    atomIsoSASAOutput$data = atomIsoSASAOutputData()
  })

  ## o5.2.1 residue SASA
  residue_sasa_null.df = data.frame(
                          ResidNe = character(),
                          Chain = character(),
                          ResidNr = integer(),
                          iCode = character(),
                          Phob.A.2 = double(),
                          Phil.A.2 = double(),
                          SASA.A.2 = double(),
                          Q.SASA. = double(),
                          N.overl. = integer(),
                          Surf.A.2 = double()
  )
  residueSASAOutput = reactiveValues(highlight = NULL, data = NULL)
  residueSASAOutputData = reactiveFileReader(2000, session,
                         function() file.path(runDir(), "id.rpopsResidue"),
                         readSasaTable, empty = residue_sasa_null.df)
  output$popsSASAResidue = DT::renderDataTable({
    residueSASAOutput$data = residueSASAOutputData()
  })

  ## o5.2.2 residue DeltaSASA
  residue_deltasasa_null.df = data.frame(
                                ResidNe = character(),
                                Chain = character(),
                                ResidNr = integer(),
                                iCode = character(),
                                D_Phob.A.2 = double(),
                                D_Phil.A.2 = double(),
                                D_SASA.A.2 = double()
  )
  residueDeltaSASAOutput = reactiveValues(highlight = NULL, data = NULL)
  residueDeltaSASAOutputData = reactiveFileReader(2000, session,
                         function() file.path(runDir(), "deltaSASA.rpopsResidue"),
                         readSasaTable, empty = residue_deltasasa_null.df)
  output$popsDeltaSASAResidue = DT::renderDataTable({
    residueDeltaSASAOutput$data = residueDeltaSASAOutputData()
  })

  ## o5.2.3 residue isolated-chain SASA
  residue_isosasa_null.df = data.frame(
                              ResidNe = character(),
                              Chain = character(),
                              ResidNr = integer(),
                              iCode = character(),
                              Phob.A.2 = double(),
                              Phil.A.2 = double(),
                              SASA.A.2 = double(),
                              Q.SASA. = double(),
                              N.overl. = integer(),
                              Surf.A.2 = double()
  )
  residueIsoSASAOutput = reactiveValues(highlight = NULL, data = NULL)
  residueIsoSASAOutputData = reactiveFileReader(2000, session,
                         function() file.path(runDir(), "isoSASA.rpopsResidue"),
                         readSasaTable, empty = residue_isosasa_null.df)
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
  chainSASAOutput = reactiveValues(highlight = NULL, data = NULL)
  chainSASAOutputData = reactiveFileReader(2000, session,
                         function() file.path(runDir(), "id.rpopsChain"),
                         readSasaTable, empty = chain_sasa_null.df)
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
  chainDeltaSASAOutput = reactiveValues(highlight = NULL, data = NULL)
  chainDeltaSASAOutputData = reactiveFileReader(2000, session,
                         function() file.path(runDir(), "deltaSASA.rpopsChain"),
                         readSasaTable, empty = chain_deltasasa_null.df)
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
  chainIsoSASAOutput = reactiveValues(highlight = NULL, data = NULL)
  chainIsoSASAOutputData = reactiveFileReader(2000, session,
                         function() file.path(runDir(), "isoSASA.rpopsChain"),
                         readSasaTable, empty = chain_isosasa_null.df)
  output$popsIsoSASAChain = DT::renderDataTable({
    chainIsoSASAOutput$data = chainIsoSASAOutputData()
  })

  ## o5.4.1 molecule SASA
  molecule_sasa_null.df = data.frame(
                            Phob.A.2 = double(),
                            Phil.A.2 = double(),
                            SASA.A.2 = double()
  )
  moleculeSASAOutput = reactiveValues(highlight = NULL, data = NULL)
  moleculeSASAOutputData = reactiveFileReader(2000, session,
                         function() file.path(runDir(), "id.rpopsMolecule"),
                         readSasaTable, empty = molecule_sasa_null.df)
  output$popsSASAMolecule = DT::renderDataTable({
    moleculeSASAOutput$data = moleculeSASAOutputData()
  })

  ## o5.4.2 molecule DeltaSASA
  molecule_deltasasa_null.df = data.frame(
                                D_Phob.A.2 = double(),
                                D_Phil.A.2 = double(),
                                D_SASA.A.2 = double()
  )
  moleculeDeltaSASAOutput = reactiveValues(highlight = NULL, data = NULL)
  moleculeDeltaSASAOutputData = reactiveFileReader(2000, session,
                         function() file.path(runDir(), "deltaSASA.rpopsMolecule"),
                         readSasaTable, empty = molecule_deltasasa_null.df)
  output$popsDeltaSASAMolecule = DT::renderDataTable({
    moleculeDeltaSASAOutput$data = moleculeDeltaSASAOutputData()
  })

  ## d1 function removed
  ## d2 download all results
  output$downloadAllResults <- downloadHandler(
    filename = function() {
      paste0("POPScomp_", if (is.null(runid_done())) "results" else runid_done(), ".zip")
    },
    content = function(file) {
      ## the identifier of the completed run, not the one of the current input
      ##   settings: those change as soon as a form field is touched
      validate(need(! is.null(runid_done()), message = "No results to download yet"))
      zipFile = file.path(sessionDir, paste0("POPScomp_", runid_done(), ".zip"))
      validate(need(file.exists(zipFile), message = "No results to download yet"))
      if (! file.copy(zipFile, file, overwrite = TRUE)) {
        stop("Could not read the result archive")
      }
    },
    contentType = "application/zip"
  )

  ## d3.1 download atom SASA
  output$downloadAtomSASA <- downloadHandler(
    filename = function() {
      paste('atomSASA_', if (is.null(runid_done())) 'results' else runid_done(), '.csv', sep = '')
    },
    content = function(fname) {
      write.csv(atomSASAOutputData(), fname)
    }
  )

  ## d3.2 download atom DeltaSASA
  output$downloadAtomDeltaSASA <- downloadHandler(
    filename = function() {
      paste('atomDeltaSASA_', if (is.null(runid_done())) 'results' else runid_done(), '.csv', sep = '')
    },
    content = function(fname) {
      write.csv(atomDeltaSASAOutputData(), fname)
    }
  )

  ## d3.3 download atom isolated-chains SASA
  output$downloadAtomIsoSASA <- downloadHandler(
    filename = function() {
      paste('atomIsoSASA_', if (is.null(runid_done())) 'results' else runid_done(), '.csv', sep = '')
    },
    content = function(fname) {
      write.csv(atomIsoSASAOutputData(), fname)
    }
  )

  ## d4.1 download residue SASA
  output$downloadResidueSASA <- downloadHandler(
    filename = function() {
      paste('residueSASA_', if (is.null(runid_done())) 'results' else runid_done(), '.csv', sep = '')
    },
    content = function(fname) {
      write.csv(residueSASAOutputData(), fname)
    }
  )

  ## d4.2 download residue DeltaSASA
  output$downloadResidueDeltaSASA <- downloadHandler(
    filename = function() {
      paste('residueDeltaSASA_', if (is.null(runid_done())) 'results' else runid_done(), '.csv', sep = '')
    },
    content = function(fname) {
      write.csv(residueDeltaSASAOutputData(), fname)
    }
  )

  ## d4.1 download residue isolated-chains SASA
  output$downloadResidueIsoSASA <- downloadHandler(
    filename = function() {
      paste('residueIsoSASA_', if (is.null(runid_done())) 'results' else runid_done(), '.csv', sep = '')
    },
    content = function(fname) {
      write.csv(residueIsoSASAOutputData(), fname)
    }
  )

  ## d5.1 download chain SASA
  output$downloadChainSASA <- downloadHandler(
    filename = function() {
      paste('chainSASA_', if (is.null(runid_done())) 'results' else runid_done(), '.csv', sep = '')
    },
    content = function(fname) {
      write.csv(chainSASAOutputData(), fname)
    }
  )

  ## d5.2 download chain DeltaSASA
  output$downloadChainDeltaSASA <- downloadHandler(
    filename = function() {
      paste('chainDeltaSASA_', if (is.null(runid_done())) 'results' else runid_done(), '.csv', sep = '')
    },
    content = function(fname) {
      write.csv(chainDeltaSASAOutputData(), fname)
    }
  )

  ## d5.1 download chain isolated-chains SASA
  output$downloadChainIsoSASA <- downloadHandler(
    filename = function() {
      paste('chainIsoSASA_', if (is.null(runid_done())) 'results' else runid_done(), '.csv', sep = '')
    },
    content = function(fname) {
      write.csv(chainIsoSASAOutputData(), fname)
    }
  )

  ## d6 download molecule SASA
  output$downloadMoleculeSASA <- downloadHandler(
    filename = function() {
      paste('moleculeSASA_', if (is.null(runid_done())) 'results' else runid_done(), '.csv', sep = '')
    },
    content = function(fname) {
      write.csv(moleculeSASAOutputData(), fname)
    }
  )

  ## d7 download molecule DeltaSASA
  output$downloadMoleculeDeltaSASA <- downloadHandler(
    filename = function() {
      paste('moleculeDeltaSASA_', if (is.null(runid_done())) 'results' else runid_done(), '.csv', sep = '')
    },
    content = function(fname) {
      write.csv(moleculeDeltaSASAOutputData(), fname)
    }
  )

  ## Readme panel text output (loaded as variable 'readme' at top of App)
  #output$readme = renderText(readme)
  #output$readme = renderUI({
  #  HTML(markdown::markdownToHTML(
  #    text = readme,
  #    fragment.only = TRUE
  #  ))
  #})
}

#_______________________________________________________________________________
# run the Shiny app
shinyApp(ui, server)

#===============================================================================
