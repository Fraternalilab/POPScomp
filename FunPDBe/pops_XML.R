#! /usr/bin/R
#===============================================================================
# run POPS ovar all XML files of PDB databank
# (C) 2018 Jens Kleinjung
#===============================================================================

library("rslurm");
#library("parallel");

#_______________________________________________________________________________
#' sadata: An S4 class input/output control.
#' @slot dirnames : array of directory names
#' @slot filenames : list of input files, list elements are the directories
#' @slot filenames_m : matrix of directory and file names in rows
#' @slot outpath : list of output paths
#' @slot command : command lines for shell
ioctrl <- setClass(
  "ioctrl",
  
  slots = c(
    dirnames = "character",
    filenames = "list",
    filenames_m = "matrix",
    outpath = "character",
    command = "character"
  )
)

ioctrl_o = ioctrl();

#_______________________________________________________________________________
## initialise resources
## synchronise XML database
#system("rsync -rlpt -v -z --delete rsync.ebi.ac.uk::pub/databases/pdb/data/structures/divided/XML/ ./XML");

## link POPS and database
#system("ln -s ~/POPS/src/pops");
#system("ln -s ~/database/XML/");

#_______________________________________________________________________________
## configure runs
## get input names of all subdirectories and 'xml.gz' files
ioctrl_o@dirnames = list.dirs(path = "./XML", full.names = FALSE);
ioctrl_o@filenames = lapply(ioctrl_o@dirnames, function(x) {
			list.files(path = paste("./XML/", x, sep = ""),
			          full.names = FALSE, pattern = 'xml\\.gz$');
})
names(ioctrl_o@filenames) = ioctrl_o@dirnames;

#_______________________________________________________________________________
## function to coerce list of directories and filenames to matrix
coerce_filenames = function(ioctrl_o) {
  dir_file_l = sapply(names(ioctrl_o@filenames), function(x) {
    if (! identical(x, "")) {
      sapply(1:length(ioctrl_o@filenames[[x]]), function(y) {
        dir_file = c(x, ioctrl_o@filenames[[x]][y]);
        return(dir_file);
      });
    }
  })
  dir_file_m = unlist(dir_file_l);
  dim(dir_file_m) = c(2, (length(dir_file_m) / 2));
  rownames(dir_file_m) = c("dir", "file");
  return(dir_file_m);
}

## populate matrix of directory and file names
ioctrl_o@filenames_m = coerce_filenames(ioctrl_o);

#_______________________________________________________________________________
## create output directory structure
## each input file will have its own output directory to accommodate multiple
##   output files from POPS
#ioctrl_o@outpath = apply(ioctrl_o@filenames_m, 2, function(x) {
#  outpath = paste("./JSON", x[1], x[2], sep = "/");
#  dir.create(paste(outpath, showWarnings = FALSE, recursive = TRUE));
#  return(outpath);
#})

#_______________________________________________________________________________
## create POPS commands
ioctrl_o@command = apply(ioctrl_o@filenames_m, 2, function(x) {
  infile = paste("./XML", x[1], x[2], sep = "/");
  outdir = paste("./JSON", x[1], x[2], sep = "/");
  ## shell command for POPSing current input file
  command = paste("./pops --pdbml", infile, "--outDirName", outdir, "--zipped --jsonOut --silent || exit 1"); 
  return(command);
})

#_______________________________________________________________________________
#_______________________________________________________________________________
## Option 1: run all command lines in serial mode
run_results = parSapply(clu, 1:length(ioctrl_o@filenames_m), function(x) {
  try(system(ioctrl_o@command[x]));
})

#_______________________________________________________________________________
## Option 2: run all command lines on Slurm cluster
## see 'rslurm' vignette: 'Parallelize R code on a Slurm cluster'

#_______________________________________________________________________________
## Option 3: run all command lines using parallelism
## number of cores
n_cores = detectCores() - 1;
## initiate cluster
clu = makeCluster(n_cores);
clusterExport(clu, "ioctrl_o");

run_results = parSapply(clu, 1:length(ioctrl_o@filenames_m), function(x) {
  try(system(ioctrl_o@command[x]));
})

stopCluster(clu);

#_______________________________________________________________________________
## validate all

#_______________________________________________________________________________
## upload all

#===============================================================================


