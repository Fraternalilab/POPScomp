
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
- A macro creates a list of all structure files in the current
    (sub)directory and parses *.xml.gz into *.json dependencies.
- POPS is invoked for each dependency file.
- POPS errors are ignored because of the leading '-' in the command call.

## Parallelism
Make is run in parallel by using the '-j' flag:
```
make -j 8 mkpops 2>&1 | tee mkpops.`date +%Y%m%d.%H%M`.log
```

