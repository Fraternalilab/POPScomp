#! /bin/bash
#===============================================================================
# Validate JSON output against funpdbe server
# see also ~/POPS/funpdbe-client/README.md
#===============================================================================

../funpdbe-client/funpdbe_client.py --path=$1 --mode=validate

