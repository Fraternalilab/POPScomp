#! /bin/sh

echo "--------------------------------------------------------------"
echo " test2a                                                       "
echo "--------------------------------------------------------------"

../src/pops --pdb 1f3r.pdb --compositionOut --typeOut --topologyOut --atomOut --residueOut --chainOut || exit 1

