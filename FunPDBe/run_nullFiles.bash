#! /bin/bash
# log the output of all nullFiles
cat nullFiles_200106.dat | xargs -i bash -c "pops --pdbml /mnt/databases/XML/{} 2>&1 | tee /mnt/databases/JSON/{}.log"

