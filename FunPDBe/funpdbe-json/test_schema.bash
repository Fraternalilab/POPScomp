#! /bin/bash
## validate example schema
python3 ./funpdbe_client.py --path=../funpdbe-json/funpdbe_example.json --mode=validate
## validate POPS atom schema
python3 ./funpdbe_client.py --path=../data/pops.pdb.allatom.json --mode=validate
## validate POPS coarse schema
python3 ./funpdbe_client.py --path=../data/pops.pdb.coarse.json --mode=validate

