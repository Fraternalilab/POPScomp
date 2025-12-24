#! /bin/bash
docker tag popscomp:latest popscomp:v3.4
docker tag popscomp:v3.4 jkleinj150999/popscomp:newest
docker push jkleinj150999/popscomp:newest
