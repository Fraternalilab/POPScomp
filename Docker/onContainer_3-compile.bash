#! /bin/bash

R CMD INSTALL POPSR

cd POPSC ; ./bootstrap ; ./configure ; make

cd

