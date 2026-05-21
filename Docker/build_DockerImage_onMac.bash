#! /bin/bash

## on Mac
docker build --platform linux/amd64 -t popscomp .
docker run --platform linux/amd64 -p 3838:3838 popscomp 

