#! /bin/bash
#===============================================================================
# Upload JSON output to FunPDBe server
# see also ~/POPS/funpdbe-client/README.md
#===============================================================================

../funpdbe-client/funpdbe_client.py --user="test-popscomp" --pwd="funpdbe-test" --api=dev --mode=post --path=$1 --resource popscomp --overwrite 

