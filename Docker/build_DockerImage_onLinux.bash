#! /bin/bash

## on Linux
docker build --no-cache -t popscomp .
docker run  -p 3838:3838 popscomp

