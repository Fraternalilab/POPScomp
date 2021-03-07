#! /bin/bash
# log the output of all nullFiles
tail +36 nullFiles_200106.dat | xargs -i bash -c "pops --pdbml /mnt/databases/XMLexcluded/{} 2>&1 | tee /mnt/databases/JSON/{}.log"

