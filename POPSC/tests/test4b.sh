#! /bin/sh

echo "--------------------------------------------------------------"
echo " test4b.sh                                                    "
echo "--------------------------------------------------------------"

../src/pops --mmcif 1F3R.cif.gz --compositionOut --typeOut --topologyOut --atomOut --residueOut --chainOut --neighbourOut --interfaceOut || exit 1

