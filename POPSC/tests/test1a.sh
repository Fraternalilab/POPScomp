#! /bin/sh

echo "--------------------------------------------------------------"
echo " test1a                                                       "
echo "--------------------------------------------------------------"

../src/pops --pdb 5lff.pdb --compositionOut --typeOut --topologyOut --atomOut --residueOut --chainOut || exit 1

