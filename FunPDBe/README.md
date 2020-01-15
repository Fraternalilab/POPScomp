
# POPS all XML files of the PDB database

Run POPS over the entire PDB database (in XML format).

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
    (sub)directory and parses *.xml.gz into *.json dependencies.
- POPS is invoked for each dependency file.
- POPS errors are ignored because of the leading '-' in the command call.

## Invoked *Makefile.validate*
- Run 'make validate' using *Makefile*.
- A macro creates a list of all *JSON* files in the current
	(sub)directory) and parses \*.json into \*.jsonval dependencies.
- The python validator *funpdbe\_client.py* is invoked for each
	dependency file.
- If validation is successful, a symbolic link \*.jsonval is created
	that points to \*.json.
- Run 'make rsyncval' to create a directory (remote or local) with all
	validated source *JSON* files. This is to upload the validated
	output to FunPDBe as one data structure.

## Invoked *Makefile.upload*
- Run 'make mkupload' using *Makefile*.
- Uploads add *JSON* output via 'lftp' to FunPDBe site.

## Parallelism
Make is run in parallel by using the '-j' flag:
```
make -j 8 mkpops 2>&1 | tee mkpops.`date +%Y%m%d.%H%M`.log
```

