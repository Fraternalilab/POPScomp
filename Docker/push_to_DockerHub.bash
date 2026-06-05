#! bin/bash

#_______________________________________________________________________________
## push 'latest'-tagged image

## login to Docker
docker login

## tag image
docker tag popscomp:latest jkleinj150999/popscomp:latest

## push image
docker push jkleinj150999/popscomp:latest

## verify
docker images


#_______________________________________________________________________________
## push 'version'-tag
docker tag popscomp:latest jkleinj150999/popscomp:v3.4
docker push jkleinj150999/popscomp:v3.4

#===============================================================================

