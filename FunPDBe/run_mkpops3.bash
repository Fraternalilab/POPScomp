#! /bin/bash
  
make -j 4 mkpops3 2>&1 | tee mkpops3.`date +%Y%m%d.%H%M`.log

