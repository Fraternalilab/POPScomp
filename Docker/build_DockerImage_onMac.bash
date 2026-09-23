#! /bin/bash

## the C program is built from POPSC, which includes the gemmi library
## as a git submodule: without this the submodule directory is empty
## and the build fails at the compilation step
git submodule update --init || exit 1

## on Mac
docker build --platform linux/amd64 -t popscomp .
docker run --platform linux/amd64 -p 3838:3838 popscomp 

