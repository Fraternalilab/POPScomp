
# POPS all XML files of the PDB database

Run POPS over the entire PDB database (in XML format).
The script *pops_XML.R* can be used to initialise, configure and run
the jobs:
- The *pops* program is symbolically linked.
- The *XML* version of the PDB database is being downloaded.
- The *JSON* output directory structure is being created such that
  each file has its own directory. Generally, this arrangement of data
  has proven useful in previous data pipelines, particularly when
  several output fileas are created per input file.
- All *xml* files are popsed via a system command.
 
