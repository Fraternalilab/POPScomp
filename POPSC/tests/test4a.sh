#! /bin/sh

echo "--------------------------------------------------------------"
echo " test4a.sh                                                    "
echo "--------------------------------------------------------------"

../src/pops --mmcif 1F3R.cif --compositionOut --typeOut --topologyOut --atomOut --residueOut --chainOut --neighbourOut --interfaceOut || exit 1

