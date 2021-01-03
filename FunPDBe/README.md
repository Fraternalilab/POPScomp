
# POPS all XML files of the PDB database
Run POPS over the entire PDB database (in XML format).

## Weekly run
The pipeline is run weekly via a cron job, because the PDB database
is updated weekly, on the virtual machine 'ace' hosted by the Synology computer 'FFJK'.
Some larger PDB files create a stack buffer overflow on 'FFJK', which is not
reproducible on a laptop running the same software. Also 'valgrind' does not
show any memory leaks. For this type of PDB structure, occasionally the pipeline
needs to be manually run on a laptop and the access mode of the created JSON files
should be adjusted via the 'make reown' step of the pipeline.

## Directory layout
- PDBML : clone of PDB database in *XML* format containing \*.xml.gz
- JSON : contains POPS output in *JSON* format
- JSONVAL : contains symbolic links to validated *JSON* output files

These directories live under *DBDIR*, which is different from
the program path *FUNPDBEDIR*, both of which can be set in the
master Makefile. 

For *make* it would be easier to have all files, i.e. PDB input,
*JSON* and symbolic link in the same directory, but for running weekly
PDB updates (rsync --delete), *lftp* uploads and general debugging purposes
the layout in separate directories seems to work better.

## Master *Makefile*
- The *XML* version of the PDB database (PDBML) is updated.
- The *JSON* output directory structure is being created 
    and the subdirectory structure of *PDBML* is mirrored.
- All structure files in PDBML are symbolically linked in
    the *JSON* directory.
- Loops over all subdirectories, invoking the *Makefile.subdir*.

## Invoked *Makefile.subdir*
- Run various targets using *Makefile*.
- A macro creates a list of all structure files in the current
    (sub)directory and parses \*.xml.gz into \*.json dependencies.
- POPS is invoked for each dependency file.
- POPS errors are ignored (i.e. *make* processing is not terminated)
    because of the leading '-' in the command call.
- Run quality controls: rename output files after POPS success,
    remove empty files where POPS failed and change the file access mode.

## Invoked *Makefile.validate*
- Run 'make validate' using *Makefile*.
- A macro creates a list of all *JSON* files.
- The python validator *funpdbe\_client.py* is invoked for each
    dependency file.
- If validation is successful, a symbolic link is created
    that points to the *JSON* file.
- Remove all files with a broken symbolic link.

## Invoked *Makefile.upload*
- Run 'make mkupload' using *Makefile*.
- Uploads add *JSON* output via 'lftp' to FunPDBe site.
    The 'mirror -RL' option uploads the files to which the symbolic links point. 

## Parallelism
Make is run in parallel by using the '-j' flag:
```
make -j 8 mkpops 2>&1 | tee mkpops.`date +%Y%m%d.%H%M`.log &
```
For split runs:
```
make -j 4 mkpops1 2>&1 | tee mkpops1.`date +%Y%m%d.%H%M`.log &
make -j 4 mkpops2 2>&1 | tee mkpops2.`date +%Y%m%d.%H%M`.log &
make -j 4 mkpops3 2>&1 | tee mkpops3.`date +%Y%m%d.%H%M`.log &
```

