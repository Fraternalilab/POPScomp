#! /bin/sh

echo "--------------------------------------------------------------"
echo " test3b.sh                                                    "
echo "--------------------------------------------------------------"

mpirun -n 2 -hostfile host_file ../src/pops --pdb 1aki.pdb --traj 1aki.sdtraj.gro --compositionOut --typeOut --topologyOut --atomOut --residueOut --chainOut || exit 1

