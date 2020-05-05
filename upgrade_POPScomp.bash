#! /bin/bash

## remove old pops binary
rm POPSC/src/pops

## pull new code from GitHub
git pull https://github.com/Fraternalilab/POPScomp

## compile POPSC
cd POPSC ; ./bootstrap ; ./configure ; make ; cd ..

## compile POPSR
R CMD build POPSR
ls -1 POPSR_*.tar.gz | tail -1 | xargs -i bash -c "R CMD INSTALL {}"

