#! /bin/bash
  
make -j 4 mkpops1 2>&1 | tee mkpops1.`date +%Y%m%d.%H%M`.log

