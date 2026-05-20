#! /bin/sh

echo "--------------------------------------------------------------"
echo " test5a.sh                                                    "
echo "--------------------------------------------------------------"

../src/pops --mmcif 5LFF.cif.gz --compositionOut --typeOut --topologyOut --atomOut --residueOut --chainOut || exit 1

