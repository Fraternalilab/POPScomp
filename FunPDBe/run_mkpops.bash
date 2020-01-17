#! /bin/bash
  
make -j 8 mkpops 2>&1 | tee mkpops.`date +%Y%m%d.%H%M`.log

