#! /bin/sh

echo "--------------------------------------------------------------"
echo " test3a                                                       "
echo "--------------------------------------------------------------"

../src/pops --pdb 1aki.pdb --traj 1aki.sdtraj.gro --compositionOut --typeOut --topologyOut --atomOut --residueOut --chainOut || exit 1

