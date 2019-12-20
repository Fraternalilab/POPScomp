#! /bin/bash
  
make -j 4 mkpops2 2>&1 | tee mkpops2.`date +%Y%m%d.%H%M`.log

