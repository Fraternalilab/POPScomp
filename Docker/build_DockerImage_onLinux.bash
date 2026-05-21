#! /bin/bash

## on Linux
docker build -t popscomp .
docker run  -p 3838:3838 popscomp

